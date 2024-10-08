import logging
import re
import zipfile
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from quantmsio.core.sdrf import SDRFHandler
from quantmsio.utils.pride_utils import get_peptidoform_proforma_version_in_mztab
from quantmsio.core.common import QUANTMSIO_VERSION, MAXQUANT_MAP, MAXQUANT_USECOLS
from quantmsio.core.feature import Feature

# format the log entries
logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)

MODIFICATION_PATTERN = re.compile(r"\((.*?\))\)")


def find_modification(peptide):
    """
    Identify the modification site based on the peptide containing modifications.

    :param peptide: Sequences of peptides
    :type peptide: str
    :return: Modification sites
    :rtype: str

    Examples:
    >>> find_modification("PEPM(UNIMOD:35)IDE")
    '4-UNIMOD:35'
    >>> find_modification("SM(UNIMOD:35)EWEIRDS(UNIMOD:21)EPTIDEK")
    '2-UNIMOD:35,9-UNIMOD:21'
    """
    peptide = str(peptide)
    original_mods = MODIFICATION_PATTERN.findall(peptide)
    peptide = MODIFICATION_PATTERN.sub(".", peptide)
    position = [i for i, x in enumerate(peptide) if x == "."]
    for j in range(1, len(position)):
        position[j] -= j

    for k in range(0, len(original_mods)):
        original_mods[k] = str(position[k]) + "-" + original_mods[k].upper()

    original_mods = ",".join(str(i) for i in original_mods) if len(original_mods) > 0 else "null"

    return original_mods


def get_mod_map(sdrf_path):
    sdrf = pd.read_csv(sdrf_path, sep="\t", nrows=1)
    mod_cols = [col for col in sdrf.columns if col.startswith("comment[modification parameters]")]
    mod_map = {}
    for col in mod_cols:
        mod_msg = sdrf[col].values[0].split(";")
        mod_dict = {k.split("=")[0]: k.split("=")[1] for k in mod_msg}
        mod = f"{mod_dict['NT']} ({mod_dict['TA']})" if "TA" in mod_dict else f"{mod_dict['NT']} ({mod_dict['PP']})"
        mod_map[mod] = mod_dict["AC"]

    return mod_map


def generate_mods(row, mod_map):
    mod_seq = row["Modified sequence"].replace("_", "")
    mod_p = find_modification(mod_seq)
    if mod_p == "null" or mod_p is None:
        return None
    for mod in row["Modifications"].split(","):
        mod = re.search(r"[A-Za-z]+.*\)$", mod)
        if mod:
            mod = mod.group()
            if mod in mod_map.keys():
                if "(" in mod_p:
                    mod_p = mod_p.replace(mod.upper(), mod_map[mod])
                else:
                    mod_p = mod_p.replace(mod[:2].upper(), mod_map[mod])
    return mod_p


class MaxQuant:
    def __init__(self, sdrf_path, evidence_path):
        self._modifications = SDRFHandler(sdrf_path).get_mods_dict()
        self._sdrf_path = sdrf_path
        self._evidence_path = evidence_path
        self.mods_map = get_mod_map(sdrf_path)

    def iter_batch(self, chunksize: int = 100000):
        for df in pd.read_csv(
            self._evidence_path,
            sep="\t",
            usecols=MAXQUANT_USECOLS + ["Potential contaminant"],
            low_memory=False,
            chunksize=chunksize,
        ):
            df = self.main_operate(df)
            yield df

    def open_from_zip_archive(self, zip_file, file_name):
        """Open file from zip archive."""
        with zipfile.ZipFile(zip_file) as z:
            with z.open(file_name) as f:
                df = pd.read_csv(f, sep="\t", usecols=MAXQUANT_USECOLS + ["Potential contaminant"], low_memory=False)
        return df

    def main_operate(self, df: pd.DataFrame):
        df.loc[:, "Modifications"] = df[["Modified sequence", "Modifications"]].apply(
            lambda row: generate_mods(row, self.mods_map), axis=1
        )
        df = df.query('`Potential contaminant`!="+"')
        df = df.drop("Potential contaminant", axis=1)
        df = df[df["PEP"] < 0.05]
        df = df.rename(columns=MAXQUANT_MAP)
        df.loc[:, "peptidoform"] = df[["sequence", "modifications"]].apply(
            lambda row: get_peptidoform_proforma_version_in_mztab(
                row["sequence"], row["modifications"], self._modifications
            ),
            axis=1,
        )
        df.loc[:, "is_decoy"] = df["is_decoy"].map({None: "0", np.nan: "0", "+": "1"})
        df["unique"] = df["pg_accessions"].apply(lambda x: "0" if ";" in str(x) else "1")
        df.loc[:, "gg_names"] = df["gg_names"].str.split(",")
        df.loc[:, "additional_scores"] = None
        df.loc[:, "modification_details"] = None
        df.loc[:, "cv_params"] = None
        df.loc[:, "quantmsio_version"] = QUANTMSIO_VERSION
        df.loc[:, "gg_accessions"] = None
        df.loc[:, "predicted_rt"] = None
        df.loc[:, "channel"] = "LFQ"
        return df

    def transform_feature(self, df: pd.DataFrame):
        df = self.merge_sdrf(df)
        df.loc[:, "global_qvalue"] = None
        df.loc[:, "pg_positions"] = None
        df.loc[:, "protein_global_qvalue"] = None
        df.loc[:, "psm_reference_file_name"] = None
        df.loc[:, "psm_scan_number"] = None
        return df

    def merge_sdrf(self, df: pd.DataFrame):
        sdrf = Feature.transform_sdrf(self._sdrf_path)
        df = pd.merge(
            df,
            sdrf,
            left_on=["reference_file_name"],
            right_on=["reference"],
            how="left",
        )
        df.drop(
            [
                "reference",
                "label",
            ],
            axis=1,
            inplace=True,
        )
        return df

    def convert_to_parquet(self, output_path: str, chunksize: int = None):
        pqwriter = None
        for df in self.iter_batch(chunksize=chunksize):
            df = self.transform_feature(df)
            feature = Feature.convert_to_parquet(df, self._modifications)
            if not pqwriter:
                pqwriter = pq.ParquetWriter(output_path, feature.schema)
            pqwriter.write_table(feature)
        if pqwriter:
            pqwriter.close()

import logging
import re
import zipfile
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from quantmsio.core.sdrf import SDRFHandler
from pyopenms import AASequence
from quantmsio.operate.tools import get_ahocorasick, get_modification_details, get_mod_map
from pyopenms.Constants import PROTON_MASS_U
from quantmsio.utils.pride_utils import get_peptidoform_proforma_version_in_mztab
from quantmsio.core.common import MAXQUANT_MAP, MAXQUANT_USECOLS
from quantmsio.core.feature import Feature
from quantmsio.core.psm import Psm

# format the log entries
logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)

MODIFICATION_PATTERN = re.compile(r"\((.*?\)?)\)")


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


def generate_mods(row, mod_map):
    mod_seq = row["Modified sequence"]
    mod_p = find_modification(mod_seq)
    if mod_p == "null" or mod_p is None:
        return None
    for mod in row["Modifications"].split(","):
        mod = re.search(r"[A-Za-z]+.*", mod)
        if mod:
            mod = mod.group()
            if mod in mod_map.keys():
                mod_p = mod_p.replace(mod.upper(), mod_map[mod])
    return mod_p


class MaxQuant:
    def __init__(self, sdrf_path, msms_path):
        self._modifications = SDRFHandler(sdrf_path).get_mods_dict()
        self._sdrf_path = sdrf_path
        self._msms_path = msms_path
        self.mods_map = get_mod_map(sdrf_path)
        self._automaton = get_ahocorasick(self.mods_map)

    def iter_batch(self, chunksize: int = 100000):
        col_df = pd.read_csv(self._msms_path, sep="\t", nrows=1)
        cols = []
        for key in self.mods_map.keys():
            col = f"{key} Probabilities"
            if col in col_df.columns:
                cols.append(col)
        for df in pd.read_csv(
            self._msms_path,
            sep="\t",
            usecols=MAXQUANT_USECOLS + cols,
            low_memory=False,
            chunksize=chunksize,
        ):
            df = self.main_operate(df)
            yield df

    def generete_calculated_mz(self, df):
        uniq_p = df["Modified sequence"].unique()
        masses_map = {k: AASequence.fromString(k).getMonoWeight() for k in uniq_p}
        mass_vector = df["Modified sequence"].map(masses_map)
        df.loc[:, "calculated_mz"] = (mass_vector + (PROTON_MASS_U * df["Charge"].values)) / df["Charge"].values

    def open_from_zip_archive(self, zip_file, file_name):
        """Open file from zip archive."""
        with zipfile.ZipFile(zip_file) as z:
            with z.open(file_name) as f:
                df = pd.read_csv(f, sep="\t", usecols=MAXQUANT_USECOLS, low_memory=False)
        return df

    def generate_modification_details(self, df):
        keys = {}
        pattern = r"\((\d+\.?\d*)\)"
        for key in self.mods_map.keys():
            col = f"{key} Probabilities"
            if col in df.columns:
                keys[key] = col
        other_mods = list(set(self.mods_map.keys()) - set(keys.keys()))

        def get_details(rows):
            modification_details = []
            for key, col in keys.items():
                modification_map = {"name": self.mods_map[key]}
                details = []
                seq = rows[col]
                if not isinstance(seq, str):
                    continue
                match_obj = re.search(pattern, seq)
                while match_obj:
                    details.append(
                        {"position": match_obj.start(0), "localization_probability": float(match_obj.group(1))}
                    )
                    seq = seq.replace(match_obj.group(0), "", 1)
                    match_obj = re.search(pattern, seq)
                modification_map["fields"] = details
                modification_details.append(modification_map)
            seq = rows["Modified sequence"]
            other_modification_details = get_modification_details(seq, self.mods_map, self._automaton, other_mods)
            modification_details = modification_details + other_modification_details
            if len(modification_details) == 0:
                return None
            else:
                return modification_details

        if len(keys.values()) != 0:
            df.loc[:, "modification_details"] = df[list(keys.values()) + ["Modified sequence"]].apply(
                get_details, axis=1
            )
        else:
            df.loc[:, "modification_details"] = None

    def main_operate(self, df: pd.DataFrame):
        df["Modified sequence"] = df["Modified sequence"].str.replace("_", "")
        df.loc[:, "Modifications"] = df[["Modified sequence", "Modifications"]].apply(
            lambda row: generate_mods(row, self.mods_map), axis=1
        )
        self.generate_modification_details(df)
        self.generete_calculated_mz(df)
        df = df[df["PEP"] < 0.05]
        df = df.rename(columns=MAXQUANT_MAP)
        df.loc[:, "peptidoform"] = df[["sequence", "modifications"]].apply(
            lambda row: get_peptidoform_proforma_version_in_mztab(
                row["sequence"], row["modifications"], self._modifications
            ),
            axis=1,
        )
        df["is_decoy"] = df["is_decoy"].map({None: "0", np.nan: "0", "+": "1"})
        df.loc[:, "unique"] = df["mp_accessions"].apply(lambda x: "0" if ";" in str(x) else "1")
        df["additional_scores"] = df["additional_scores"].apply(
            lambda x: [{"name": "maxquant", "value": np.float32(x)}]
        )
        df.loc[:, "best_id_score"] = None
        df.loc[:, "cv_params"] = None
        df.loc[:, "consensus_support"] = None
        df.loc[:, "predicted_rt"] = None
        df.loc[:, "global_qvalue"] = None
        df.loc[:, "pg_global_qvalue"] = None
        return df

    @staticmethod
    def transform_psm(df: pd.DataFrame):
        df.loc[:, "intensity_array"] = None
        df.loc[:, "ion_mobility"] = None
        df.loc[:, "mz_array"] = None
        df.loc[:, "number_peaks"] = None
        df.loc[:, "rank"] = None

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
            self.transform_psm(df)
            Psm.convert_to_parquet_format(df, self._modifications)
            parquet = Psm.transform_parquet(df)
            if not pqwriter:
                pqwriter = pq.ParquetWriter(output_path, parquet.schema)
            pqwriter.write_table(parquet)
        if pqwriter:
            pqwriter.close()

    # def convert_to_parquet(self, output_path: str, chunksize: int = None):
    #     pqwriter = None
    #     for df in self.iter_batch(chunksize=chunksize):
    #         df = self.transform_feature(df)
    #         Feature.convert_to_parquet_format(df, self._modifications)
    #         feature = Feature.transform_feature(df)
    #         if not pqwriter:
    #             pqwriter = pq.ParquetWriter(output_path, feature.schema)
    #         pqwriter.write_table(feature)
    #     if pqwriter:
    #         pqwriter.close()

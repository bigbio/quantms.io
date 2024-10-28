import logging
import re
import zipfile
import numpy as np
import pandas as pd
import codecs
import os
import pyarrow.parquet as pq
from pyopenms import ModificationsDB
from pyopenms import AASequence
from quantmsio.operate.tools import get_ahocorasick, get_modification_details
from pyopenms.Constants import PROTON_MASS_U
from quantmsio.core.common import MAXQUANT_MAP, MAXQUANT_USECOLS
from quantmsio.core.feature import Feature
from quantmsio.core.psm import Psm

logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)


class MaxQuant:
    def __init__(self, msms_path):
        self._msms_path = msms_path
        self.mods_map = self.get_mods_map()
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
            df = df.rename(columns=MAXQUANT_MAP)
            df = self.main_operate(df)
            yield df

    def generete_calculated_mz(self, df):
        uniq_p = df["peptidoform"].unique()
        masses_map = {k: AASequence.fromString(k).getMonoWeight() for k in uniq_p}
        mass_vector = df["peptidoform"].map(masses_map)
        df.loc[:, "calculated_mz"] = (mass_vector + (PROTON_MASS_U * df["precursor_charge"].values)) / df[
            "precursor_charge"
        ].values

    def open_from_zip_archive(self, zip_file, file_name):
        """Open file from zip archive."""
        with zipfile.ZipFile(zip_file) as z:
            with z.open(file_name) as f:
                df = pd.read_csv(f, sep="\t", usecols=MAXQUANT_USECOLS, low_memory=False)
        return df

    def get_mods_map(self):
        if os.stat(self._msms_path).st_size == 0:
            raise ValueError("File is empty")
        f = codecs.open(self._msms_path, "r", "utf-8")
        line = f.readline()
        pattern = r"sequence\s+(.*?)\s+Proteins"
        match = re.search(pattern, line, re.DOTALL)
        mod_seq = match.group(1)
        mods_map = {}
        modifications_db = ModificationsDB()
        if mod_seq:
            for mod in mod_seq.split("\t"):
                if "Probabilities" not in mod and "diffs" not in mod:
                    name = re.search(r"([^ ]+)\s?", mod)
                    Mod = modifications_db.getModification(name.group(1))
                    unimod = Mod.getUniModAccession()
                    match = re.search(r"\((.*?)\)", mod)
                    if match:
                        site = match.group(1)
                    else:
                        site = "X"
                    mods_map[mod] = [unimod.upper(), site]
                    mods_map[unimod.upper()] = [mod, site]
        return mods_map

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
                modification_map = {"modification_name": self.mods_map[key][0]}
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
            seq = rows["peptidoform"]
            peptidoform, other_modification_details = get_modification_details(
                seq, self.mods_map, self._automaton, other_mods
            )
            modification_details = modification_details + other_modification_details
            if len(modification_details) == 0:
                return [peptidoform, None]
            else:
                return [peptidoform, modification_details]

        df[["peptidoform", "modifications"]] = df[list(keys.values()) + ["peptidoform"]].apply(
            get_details, axis=1, result_type="expand"
        )

    def main_operate(self, df: pd.DataFrame):
        df["peptidoform"] = df["peptidoform"].str.replace("_", "")
        self.generete_calculated_mz(df)
        self.generate_modification_details(df)
        df = df[df["posterior_error_probability"] < 0.05].copy()
        df["is_decoy"] = df["is_decoy"].map({None: "0", np.nan: "0", "+": "1"})
        df["additional_scores"] = df["additional_scores"].apply(
            lambda x: [{"name": "maxquant", "value": np.float32(x)}]
        )
        # df.loc[:, "best_id_score"] = None
        df.loc[:, "cv_params"] = None
        df.loc[:, "predicted_rt"] = None
        return df

    @staticmethod
    def transform_psm(df: pd.DataFrame):
        df.loc[:, "intensity_array"] = None
        df.loc[:, "ion_mobility"] = None
        df.loc[:, "mz_array"] = None
        df.loc[:, "number_peaks"] = None

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
            Psm.convert_to_parquet_format(df)
            parquet = Psm.transform_parquet(df)
            if not pqwriter:
                pqwriter = pq.ParquetWriter(output_path, parquet.schema)
            pqwriter.write_table(parquet)
        if pqwriter:
            pqwriter.close()

import pyarrow as pa
import pyarrow.parquet as pq
from quantmsio.utils.file_utils import extract_protein_list
from quantmsio.utils.pride_utils import generate_scan_number
from quantmsio.utils.pride_utils import get_peptidoform_proforma_version_in_mztab
from quantmsio.core.common import PSM_USECOLS, PSM_MAP, QUANTMSIO_VERSION
from quantmsio.core.mztab import MzTab, generate_modification_list
from quantmsio.core.format import PSM_FIELDS
import pandas as pd

PSM_SCHEMA = pa.schema(
    PSM_FIELDS,
    metadata={"description": "psm file in quantms.io format"},
)


class Psm(MzTab):
    def __init__(self, mzTab_path):
        super(Psm, self).__init__(mzTab_path)
        self._ms_runs = self.extract_ms_runs()
        self._protein_global_qvalue_map = self.get_protein_map()
        self._modifications = self.get_modifications()
        self._score_names = self.get_score_names()

    def iter_psm_table(self, chunksize=1000000, protein_str=None):
        for df in self.skip_and_load_csv("PSH", chunksize=chunksize):
            if protein_str:
                df = df[df["accession"].str.contains(f"{protein_str}", na=False)]
            no_cols = set(PSM_USECOLS) - set(df.columns)
            for col in no_cols:
                if col == "unique":
                    df.loc[:, col] = df["accession"].apply(lambda x: 0 if ";" in x else 1)
                else:
                    df.loc[:, col] = None
            df.rename(columns=PSM_MAP, inplace=True)
            yield df

    def generate_report(self, chunksize=1000000, protein_str=None):
        for df in self.iter_psm_table(chunksize=chunksize, protein_str=protein_str):
            self.transform_psm(df)
            self.add_addition_msg(df)
            df = self.convert_to_parquet(df, self._modifications)
            yield df

    def transform_psm(self, df):
        df.loc[:, "pg_positions"] = df[["start", "end"]].apply(
            lambda row: self.generate_positions(row["start"], row["end"]), axis=1
        )
        df.loc[:, "scan_number"] = df["spectra_ref"].apply(lambda x: generate_scan_number(x))

        df.loc[:, "reference_file_name"] = df["spectra_ref"].apply(lambda x: self._ms_runs[x[: x.index(":")]])
        df.loc[:, "additional_scores"] = df[list(self._score_names.values())].apply(
            self._genarate_additional_scores, axis=1
        )
        df.loc[:, "peptidoform"] = df[["modifications", "sequence"]].apply(
            lambda row: get_peptidoform_proforma_version_in_mztab(
                row["sequence"], row["modifications"], self._modifications
            ),
            axis=1,
        )
        df.drop(["start", "end", "spectra_ref", "search_engine", "search_engine_score[1]"], inplace=True, axis=1)

    def _genarate_additional_scores(self, cols):
        struct_list = []
        for software, score in self._score_names.items():
            struct = {"name": software, "value": cols[score]}
            struct_list.append(struct)
        return struct_list

    def add_addition_msg(self, df):
        df.loc[:, "protein_global_qvalue"] = df["pg_accessions"].map(self._protein_global_qvalue_map)
        df.loc[:, "modification_details"] = None
        df.loc[:, "predicted_rt"] = None
        df.loc[:, "gg_accessions"] = None
        df.loc[:, "gg_names"] = None
        df.loc[:, "ion_mobility"] = None
        df.loc[:, "num_peaks"] = None
        df.loc[:, "mz_array"] = None
        df.loc[:, "intensity_array"] = None
        df.loc[:, "rank"] = None
        df.loc[:, "cv_params"] = None
        df.loc[:, "quantmsio_version"] = QUANTMSIO_VERSION

    def write_feature_to_file(self, output_path, chunksize=1000000, protein_file=None):
        protein_list = extract_protein_list(protein_file) if protein_file else None
        protein_str = "|".join(protein_list) if protein_list else None
        pqwriter = None
        for p in self.generate_report(chunksize=chunksize, protein_str=protein_str):
            if not pqwriter:
                pqwriter = pq.ParquetWriter(output_path, p.schema)
            pqwriter.write_table(p)
        if pqwriter:
            pqwriter.close()

    @staticmethod
    def convert_to_parquet(res, modifications):
        res["pg_accessions"] = res["pg_accessions"].str.split(";")
        res["protein_global_qvalue"] = res["protein_global_qvalue"].astype(float)
        res["unique"] = res["unique"].astype("Int32")
        res["modifications"] = res["modifications"].apply(lambda x: generate_modification_list(x, modifications))
        res["precursor_charge"] = res["precursor_charge"].map(lambda x: None if pd.isna(x) else int(x)).astype("Int32")
        res["calculated_mz"] = res["calculated_mz"].astype(float)
        res["observed_mz"] = res["observed_mz"].astype(float)
        res["posterior_error_probability"] = res["posterior_error_probability"].astype(float)
        res["global_qvalue"] = res["global_qvalue"].astype(float)
        res["is_decoy"] = res["is_decoy"].map(lambda x: None if pd.isna(x) else int(x)).astype("Int32")

        res["scan_number"] = res["scan_number"].astype(str)

        if "rt" in res.columns:
            res["rt"] = res["rt"].astype(float)
        else:
            res.loc[:, "rt"] = None
        return pa.Table.from_pandas(res, schema=PSM_SCHEMA)

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from quantmsio.core.feature_in_memory import FeatureInMemory
from quantmsio.core.feature_in_memory import get_modifications
from quantmsio.utils.pride_utils import generate_scan_number
from quantmsio.utils.pride_utils import get_peptidoform_proforma_version_in_mztab
from quantmsio.utils.pride_utils import get_petidoform_msstats_notation


class PsmInMemory:
    def __init__(self, schema):
        self.schema = schema
        self._psms_columns = None
        self._ms_runs = None
        self._score_names = None
        self._feature = FeatureInMemory()
        self._psm_usecols = [
            "sequence",
            "accession",
            "unique",
            "modifications",
            "retention_time",
            "charge",
            "exp_mass_to_charge",
            "calc_mass_to_charge",
            "spectra_ref",
            "start",
            "end",
            "opt_global_cv_MS:1000889_peptidoform_sequence",
            "opt_global_cv_MS:1002217_decoy_peptide",
            "opt_global_Posterior_Error_Probability_score",
            "opt_global_consensus_support",
            "search_engine_score[1]",
        ]
        self._psm_map = {
            "accession": "protein_accessions",
            "start": "protein_start_positions",
            "end": "protein_end_positions",
            "opt_global_cv_MS:1000889_peptidoform_sequence": "peptidoform",
            "opt_global_Posterior_Error_Probability_score": "posterior_error_probability",
            "opt_global_cv_MS:1002217_decoy_peptide": "is_decoy",
            "spectra_ref": "reference_file_name",
            "opt_global_consensus_support": "consensus_support",
        }

    def generate_psm_parquet(self, mztab_path, chunksize=1000000):
        protein_map = self._feature._get_protein_map(mztab_path)
        self._ms_runs = self._feature.extract_ms_runs(mztab_path)
        self._feature._ms_runs = self._ms_runs
        self._score_names = self._feature._get_score_names(mztab_path)
        map_dict = self._feature._extract_from_pep(mztab_path)
        psms = self._feature.skip_and_load_csv(
            mztab_path,
            "PSH",
            sep="\t",
            dtype={"start": str, "end": str},
            chunksize=chunksize,
        )
        self._psms_columns = self._feature._psms_columns
        self._modifications = get_modifications(mztab_path)
        for psm in psms:
            if "opt_global_cv_MS:1000889_peptidoform_sequence" not in psm.columns:
                psm.loc[:, "opt_global_cv_MS:1000889_peptidoform_sequence"] = psm[["modifications", "sequence"]].apply(
                    lambda row: get_petidoform_msstats_notation(
                        row["sequence"], row["modifications"], self._modifications
                    ),
                    axis=1,
                )
            if "opt_global_cv_MS:1002217_decoy_peptide" not in psm.columns:
                psm.loc[:, "opt_global_cv_MS:1002217_decoy_peptide"] = "0"
            for col in self._psm_usecols:
                if col not in psm.columns:
                    psm.loc[:, col] = None
            psm.loc[:, "scan_number"] = psm["spectra_ref"].apply(lambda x: generate_scan_number(x))
            psm["spectra_ref"] = psm["spectra_ref"].apply(lambda x: self._ms_runs[x.split(":")[0]])
            self._psm_usecols.append("scan_number")
            psm = psm[self._psm_usecols]
            self._psm_usecols.pop()
            psm.rename(columns=self._psm_map, inplace=True)
            psm.loc[:, "global_qvalue"] = psm[["peptidoform", "charge"]].apply(
                lambda row: (
                    map_dict[(row["peptidoform"], row["charge"])][0]
                    if (row["peptidoform"], row["charge"]) in map_dict
                    else None
                ),
                axis=1,
            )
            psm.loc[:, "protein_global_qvalue"] = psm["protein_accessions"].apply(
                lambda x: self._feature._handle_protein_map(protein_map, x)
            )
            psm_score_name = self._score_names["psm_score"]
            psm.loc[:, "id_scores"] = (
                psm_score_name
                + ":"
                + psm["search_engine_score[1]"].astype(str)
                + ","
                + "Best PSM PEP:"
                + psm["posterior_error_probability"].astype(str)
            )
            psm.drop(["search_engine_score[1]"], inplace=True, axis=1)
            psm.loc[:, "peptidoform"] = psm[["sequence", "modifications"]].apply(
                lambda row: get_peptidoform_proforma_version_in_mztab(
                    row["sequence"], row["modifications"], self._modifications
                ),
                axis=1,
            )
            parquet_table = self.convert_to_parquet(psm)
            yield parquet_table

    def write_feature_to_file(self, mztab_path, output_path, chunksize=1000000):
        """
        write parquet to file
        """
        pqwriter = None
        for feature in self.generate_psm_parquet(mztab_path, chunksize=chunksize):
            if not pqwriter:
                pqwriter = pq.ParquetWriter(output_path, feature.schema)
            pqwriter.write_table(feature)
        if pqwriter:
            pqwriter.close()

    @staticmethod
    def __split_start_or_end(value):
        """
        split start or end
        :param value: start or end
        """
        if pd.isna(value) or value == "null":
            return pd.NA
        elif "," in str(value):
            return list(map(int, value.split(",")))
        elif value is np.nan:
            return None
        else:
            return [int(value)]

    def convert_to_parquet(self, res):
        """
        res: psm dataframe
        return: parquet table
        """
        self._feature._modifications = self._modifications
        self._feature._score_names = self._score_names
        res["id_scores"] = res["id_scores"].apply(lambda x: x.split(","))
        res["sequence"] = res["sequence"].astype(str)
        res["protein_accessions"] = res["protein_accessions"].str.split(";")
        res["protein_start_positions"] = res["protein_start_positions"].apply(self.__split_start_or_end).to_list()
        res["protein_end_positions"] = res["protein_end_positions"].apply(self.__split_start_or_end).to_list()
        res["protein_global_qvalue"] = res["protein_global_qvalue"].astype(float)
        res["unique"] = res["unique"].map(lambda x: pd.NA if pd.isna(x) else int(x)).astype("Int32")
        res["modifications"] = res["modifications"].apply(lambda x: self._feature._generate_modification_list(x))
        res["charge"] = res["charge"].map(lambda x: pd.NA if pd.isna(x) else int(x)).astype("Int32")
        res["exp_mass_to_charge"] = res["exp_mass_to_charge"].astype(float)
        res["calc_mass_to_charge"] = res["calc_mass_to_charge"].astype(float)
        res["posterior_error_probability"] = res["posterior_error_probability"].astype(float)
        res["global_qvalue"] = res["global_qvalue"].astype(float)
        res["is_decoy"] = res["is_decoy"].map(lambda x: pd.NA if pd.isna(x) else int(x)).astype("Int32")
        res["consensus_support"] = res["consensus_support"].astype(float)
        res["scan_number"] = res["scan_number"].astype(str)

        if "retention_time" in res.columns:
            res["retention_time"] = res["retention_time"].astype(float)
        else:
            res.loc[:, "retention_time"] = None

        res.loc[:, "num_peaks"] = None
        res.loc[:, "mz_array"] = None
        res.loc[:, "intensity_array"] = None

        res.loc[:, "gene_accessions"] = None
        res.loc[:, "gene_names"] = None
        return pa.Table.from_pandas(res, schema=self.schema)

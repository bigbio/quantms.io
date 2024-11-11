import time
import logging
import numpy as np
import pandas as pd
import os
import pyarrow.parquet as pq
import concurrent.futures
from pathlib import Path
from pyopenms import AASequence
from pyopenms.Constants import PROTON_MASS_U
from quantmsio.operate.tools import get_ahocorasick
from quantmsio.utils.file_utils import close_file, extract_protein_list, save_slice_file
from quantmsio.core.sdrf import SDRFHandler
from quantmsio.core.mztab import MzTab
from quantmsio.core.feature import Feature
from quantmsio.core.duckdb import DuckDB
from quantmsio.utils.pride_utils import generate_scan_number
from quantmsio.core.common import DIANN_MAP


class DiaNNConvert(DuckDB):

    def __init__(self, diann_report, sdrf_path, duckdb_max_memory="16GB", duckdb_threads=4):
        super(DiaNNConvert, self).__init__(diann_report, duckdb_max_memory, duckdb_threads)
        self.init_duckdb()
        self._sdrf = SDRFHandler(sdrf_path)
        self._mods_map = self._sdrf.get_mods_dict()
        self._automaton = get_ahocorasick(self._mods_map)
        self._sample_map = self._sdrf.get_sample_map_run()

    def init_duckdb(self):

        self._duckdb.execute("""CREATE INDEX idx_precursor_q ON report ("Precursor.Id", "Q.Value")""")
        self._duckdb.execute("""CREATE INDEX idx_run ON report ("Run")""")

    def get_report_from_database(self, runs: list) -> pd.DataFrame:
        """
        This function loads the report from the duckdb database for a group of ms_runs.
        :param runs: A list of ms_runs
        :return: The report
        """
        s = time.time()
        database = self._duckdb.query(
            """
            select "File.Name", "Run", "Protein.Group", "Protein.Ids", "Genes", "RT", "Predicted.RT", "RT.Start", "RT.Stop", "Precursor.Id", "Q.Value", "Global.Q.Value", "PEP",
                   "PG.Q.Value", "Global.PG.Q.Value", "Modified.Sequence", "Stripped.Sequence", "Precursor.Charge", "Precursor.Quantity", "Precursor.Normalised"
                    from report
            where Run IN {}
            """.format(
                tuple(runs)
            )
        )
        report = database.df()
        et = time.time() - s
        logging.info("Time to load report {} seconds".format(et))
        return report

    def get_masses_and_modifications_map(self):
        database = self._duckdb.query(
            """
            select DISTINCT "Modified.Sequence" from report
            """
        )
        report = database.df()
        uniq_p = report["Modified.Sequence"].values
        masses_map = {k: AASequence.fromString(k).getMonoWeight() for k in uniq_p}
        modifications_map = {k: AASequence.fromString(k).toString() for k in uniq_p}

        return masses_map, modifications_map

    def get_peptide_map_from_database(self):
        s = time.time()
        database = self._duckdb.query(
            """
            SELECT "Precursor.Id","Q.Value","Run"
            FROM (
            SELECT "Precursor.Id", "Q.Value","Run", ROW_NUMBER()
            OVER (PARTITION BY "Precursor.Id" ORDER BY "Q.Value" ASC) AS row_num
            FROM report
            ) AS subquery
            WHERE row_num = 1;
            """
        )
        peptide_df = database.df()
        peptide_df.set_index("Precursor.Id", inplace=True)
        # peptide_map = peptide_df.to_dict()["Q.Value"]
        best_ref_map = peptide_df.to_dict()["Run"]
        et = time.time() - s
        logging.info("Time to load peptide map {} seconds".format(et))
        return best_ref_map

    def main_report_df(self, qvalue_threshold: float, mzml_info_folder: str, file_num: int, protein_str: str = None):
        def intergrate_msg(n):
            nonlocal report
            nonlocal mzml_info_folder
            files = list(Path(mzml_info_folder).glob(f"*{n}_mzml_info.tsv"))
            if not files:
                raise ValueError(f"Could not find {n} info file in {dir}")
            target = pd.read_csv(
                files[0],
                sep="\t",
                usecols=["Retention_Time", "SpectrumID", "Exp_Mass_To_Charge"],
            )
            group = report[report["Run"] == n].copy()
            group.sort_values(by="rt_start", inplace=True)
            target.rename(
                columns={
                    "Retention_Time": "rt_start",
                    "SpectrumID": "scan",
                    "Exp_Mass_To_Charge": "observed_mz",
                },
                inplace=True,
            )
            target["rt_start"] = target["rt_start"] / 60
            res = pd.merge_asof(group, target, on="rt_start", direction="nearest")
            return res

        # query duckdb
        # best_ref_map = self.get_peptide_map_from_database()
        masses_map, modifications_map = self.get_masses_and_modifications_map()

        info_list = [
            mzml.replace("_mzml_info.tsv", "")
            for mzml in os.listdir(mzml_info_folder)
            if mzml.endswith("_mzml_info.tsv")
        ]
        info_list = [info_list[i : i + file_num] for i in range(0, len(info_list), file_num)]
        for refs in info_list:
            report = self.get_report_from_database(refs)
            report.rename(columns=DIANN_MAP, inplace=True)
            if protein_str:
                report = report[report["pg_accessions"].str.contains(f"{protein_str}", na=False)]
            # restrict
            report = report[report["qvalue"] < qvalue_threshold]
            usecols = report.columns.to_list() + [
                "scan",
                "observed_mz",
            ]
            with concurrent.futures.ThreadPoolExecutor(100) as executor:
                results = executor.map(intergrate_msg, refs)
            report = np.vstack([result.values for result in results])
            report = pd.DataFrame(report, columns=usecols)

            # cal value and mod
            mass_vector = report["peptidoform"].map(masses_map)
            report["calculated_mz"] = (mass_vector + (PROTON_MASS_U * report["precursor_charge"].values)) / report[
                "precursor_charge"
            ].values

            report["peptidoform"] = report["peptidoform"].map(modifications_map)

            yield report

    def add_additional_msg(self, report: pd.DataFrame):
        """
        Perform some transformations in the report dataframe to help with the generation of the psm and feature files.
        :param report: The report dataframe
        """
        select_mods = list(self._mods_map.keys())
        report["reference_file_name"] = report["reference_file_name"].apply(lambda x: x.split(".")[0])
        report[["peptidoform", "modifications"]] = report[["peptidoform"]].apply(
            lambda row: MzTab.generate_modifications_details(
                row["peptidoform"], self._mods_map, self._automaton, select_mods
            ),
            axis=1,
            result_type="expand",
        )
        report.loc[:, "channel"] = "LFQ"
        report.loc[:, "intensities"] = report[["reference_file_name", "channel", "intensity"]].apply(
            lambda rows: [
                {
                    "sample_accession": self._sample_map[rows["reference_file_name"] + "-" + rows["channel"]],
                    "channel": rows["channel"],
                    "intensity": rows["intensity"],
                }
            ],
            axis=1,
        )
        report.loc[:, "is_decoy"] = "0"
        report.loc[:, "unique"] = report["pg_accessions"].apply(lambda x: "0" if ";" in str(x) else "1")
        report["scan"] = report["scan"].apply(generate_scan_number)
        report["mp_accessions"] = report["mp_accessions"].str.split(";")
        report["pg_accessions"] = report["pg_accessions"].str.split(";")
        report.loc[:, "anchor_protein"] = report["pg_accessions"].str[0]
        report.loc[:, "gg_names"] = report["gg_names"].str.split(",")

        report.loc[:, "additional_intensities"] = report[
            ["reference_file_name", "channel", "normalize_intensity"]
        ].apply(
            lambda rows: [
                {
                    "sample_accession": self._sample_map[rows["reference_file_name"] + "-" + rows["channel"]],
                    "channel": rows["channel"],
                    "additional_intensity": [
                        {"intensity_name": "normalize_intensity", "intensity_value": rows["normalize_intensity"]}
                    ],
                }
            ],
            axis=1,
        )
        report.loc[:, "additional_scores"] = report[["qvalue", "pg_qvalue"]].apply(
            lambda row: [
                {"score_name": "qvalue", "score_value": row["qvalue"]},
                {"score_name": "pg_qvalue", "score_value": row["pg_qvalue"]},
            ],
            axis=1,
        )
        report.loc[:, "scan_reference_file_name"] = None
        report.loc[:, "scan"] = None
        report.loc[:, "cv_params"] = None
        report.loc[:, "gg_accessions"] = None
        report.loc[:, "ion_mobility"] = None
        report.loc[:, "start_ion_mobility"] = None
        report.loc[:, "stop_ion_mobility"] = None

    def generate_feature(
        self, qvalue_threshold: float, mzml_info_folder: str, file_num: int = 50, protein_str: str = None
    ):
        for report in self.main_report_df(qvalue_threshold, mzml_info_folder, file_num, protein_str):
            s = time.time()
            self.add_additional_msg(report)
            Feature.convert_to_parquet_format(report)
            et = time.time() - s
            logging.info("Time to generate psm and feature file {} seconds".format(et))
            yield report

    def write_feature_to_file(
        self,
        qvalue_threshold: float,
        mzml_info_folder: str,
        output_path: str,
        file_num: int = 50,
        protein_file=None,
    ):
        protein_list = extract_protein_list(protein_file) if protein_file else None
        protein_str = "|".join(protein_list) if protein_list else None
        pqwriter = None
        for report in self.generate_feature(qvalue_threshold, mzml_info_folder, file_num, protein_str):
            feature = Feature.transform_feature(report)
            if not pqwriter:
                pqwriter = pq.ParquetWriter(output_path, feature.schema)
            pqwriter.write_table(feature)
        close_file(pqwriter=pqwriter)
        self.destroy_duckdb_database()

    def write_features_to_file(
        self,
        qvalue_threshold: float,
        mzml_info_folder: str,
        output_folder: str,
        filename: str,
        partitions: list,
        file_num:int = 50,
        protein_file=None,
    ):
        pqwriters = {}
        protein_list = extract_protein_list(protein_file) if protein_file else None
        protein_str = "|".join(protein_list) if protein_list else None
        for report in self.generate_feature(qvalue_threshold, mzml_info_folder, file_num, protein_str):
            for key, df in Feature.slice(report, partitions):
                feature = Feature.transform_feature(df)
                pqwriters = save_slice_file(feature, pqwriters, output_folder, key, filename)
        close_file(pqwriters=pqwriters)
        self.destroy_duckdb_database()
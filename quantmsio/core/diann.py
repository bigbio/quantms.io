import concurrent.futures
import logging
import os
import time
from collections import defaultdict
from pathlib import Path
from typing import Union

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pyopenms import AASequence
from pyopenms.Constants import PROTON_MASS_U

from quantmsio.core.common import DIANN_MAP, DIANN_USECOLS, DIANN_PG_MAP, DIANN_PG_USECOLS, PG_SCHEMA
from quantmsio.core.duckdb import DuckDB
from quantmsio.core.feature import Feature
from quantmsio.core.mztab import MzTab
from quantmsio.core.sdrf import SDRFHandler
from quantmsio.operate.tools import get_ahocorasick
from quantmsio.utils.file_utils import close_file, extract_protein_list, save_slice_file
from quantmsio.utils.pride_utils import generate_scan_number

DIANN_SQL = ", ".join([f'"{name}"' for name in DIANN_USECOLS])
DIANN_PG_SQL = ", ".join([f'"{name}"' for name in DIANN_PG_USECOLS])


class DiaNNConvert(DuckDB):

    def __init__(self, diann_report, sdrf_path=None, duckdb_max_memory="16GB", duckdb_threads=4):
        super(DiaNNConvert, self).__init__(diann_report, duckdb_max_memory, duckdb_threads)
        if sdrf_path:
            self._sdrf = SDRFHandler(sdrf_path)
            self._mods_map = self._sdrf.get_mods_dict()
            self._automaton = get_ahocorasick(self._mods_map)
            self._sample_map = self._sdrf.get_sample_map_run()

    def get_report_from_database(self, runs: list, sql: str = DIANN_SQL) -> pd.DataFrame:

        s = time.time()
        database = self._duckdb.query(
            """
            select {}
            from report
            where Run IN {}
            """.format(
                sql, tuple(runs)
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

    @staticmethod
    def get_peptide_count(df):
        peptide_count = defaultdict(int)
        for proteins in df["pg_accessions"]:
            for protein in proteins.split(";"):
                peptide_count[protein] += 1
        return peptide_count

    def generate_pg_matrix(self, report):
        peptide_count = self.get_peptide_count(report)
        report.drop_duplicates(subset=["pg_accessions"], inplace=True)
        report["gg_accessions"] = report["gg_accessions"].str.split(";")
        report["pg_names"] = report["pg_names"].str.split(";")
        report["reference_file_name"] = report["reference_file_name"].apply(lambda x: x.split(".")[0])
        report["pg_accessions"] = report["pg_accessions"].str.split(";")
        report.loc[:, "peptides"] = report["pg_accessions"].apply(
            lambda proteins: [
                {"protein_name": protein, "peptide_count": peptide_count[protein]} for protein in proteins
            ]
        )
        report.loc[:, "is_decoy"] = 0
        report.loc[:, "additional_intensities"] = report[["normalize_intensity", "lfq"]].apply(
            lambda rows: [
                {"intensity_name": name, "intensity_value": rows[name]} for name in ["normalize_intensity", "lfq"]
            ],
            axis=1,
        )
        report.loc[:, "additional_scores"] = report["qvalue"].apply(
            lambda value: [{"score_name": "qvalue", "score_value": value}]
        )
        report.loc[:, "contaminant"] = None
        report.loc[:, "anchor_protein"] = None
        return report

    def main_report_df(
        self, qvalue_threshold: float, mzml_info_folder: Union[Path, str], file_num: int, protein_str: str = None
    ):
        """
        Process DIA-NN report data and integrate with MS info.
        Optimized for better performance with large datasets.

        Parameters:
        -----------
        qvalue_threshold : float
            Threshold for filtering by q-value
        mzml_info_folder : str
            Path to folder containing MS info parquet files
        file_num : int
            Number of files to process in each batch
        protein_str : str, optional
            Protein accession filter

        Yields:
        -------
        pandas.DataFrame
            Processed report data
        """

        def integrate_msg(n):
            """Integrate MS info with report data for a specific run"""
            # Use Path for more reliable file handling
            ms_info_file = next(Path(mzml_info_folder).glob(f"*{n}_ms_info.parquet"), None)
            if not ms_info_file:
                raise ValueError(f"Could not find MS info file for run {n} in {mzml_info_folder}")

            # Read only necessary columns
            target = pd.read_parquet(ms_info_file, columns=["rt", "scan", "observed_mz"])

            # Filter report data for this run (avoid copy if possible)
            group = report_filtered[report_filtered["run"] == n]

            # Sort by retention time
            group = group.sort_values(by="rt")

            # Convert retention time to minutes
            target["rt"] = target["rt"] / 60

            # Merge using pandas merge_asof for nearest match
            res = pd.merge_asof(group, target, on="rt", direction="nearest")
            return res

        # Get mass and modification maps once
        masses_map, modifications_map = self.get_masses_and_modifications_map()

        # Get list of MS info files more efficiently
        info_files = Path(mzml_info_folder).glob("*_ms_info.parquet")
        info_list = [f.stem.replace("_ms_info", "") for f in info_files]

        # Process in batches
        for i in range(0, len(info_list), file_num):
            refs = info_list[i : i + file_num]

            # Get report data for this batch
            report = self.get_report_from_database(refs)

            # Rename columns
            report.rename(columns=DIANN_MAP, inplace=True)

            # Filter data
            report = report.dropna(subset=["pg_accessions"])
            if protein_str:
                report = report[report["pg_accessions"].str.contains(protein_str, na=False)]

            # Apply q-value threshold
            report_filtered = report[report["qvalue"] < qvalue_threshold]

            # Define columns for final dataframe
            usecols = report_filtered.columns.tolist() + ["scan", "observed_mz"]

            # Process each run in parallel with a more reasonable number of workers
            max_workers = min(len(refs), os.cpu_count() * 2)
            with concurrent.futures.ThreadPoolExecutor(max_workers) as executor:
                results = list(executor.map(integrate_msg, refs))

            # Combine results more efficiently
            if results:
                combined_report = pd.concat(results, ignore_index=True)

                # Calculate m/z values vectorized
                mass_vector = combined_report["peptidoform"].map(masses_map)
                combined_report["calculated_mz"] = (
                    mass_vector + (PROTON_MASS_U * combined_report["precursor_charge"])
                ) / combined_report["precursor_charge"]

                # Map peptidoforms
                combined_report["peptidoform"] = combined_report["peptidoform"].map(modifications_map)

                yield combined_report

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
            ["reference_file_name", "channel", "normalize_intensity", "lfq"]
        ].apply(
            lambda rows: [
                {
                    "sample_accession": self._sample_map[rows["reference_file_name"] + "-" + rows["channel"]],
                    "channel": rows["channel"],
                    "additional_intensity": [
                        {"intensity_name": "normalize_intensity", "intensity_value": rows["normalize_intensity"]},
                        {"intensity_name": "lfq", "intensity_value": rows["lfq"]},
                    ],
                }
            ],
            axis=1,
        )
        report.loc[:, "additional_scores"] = report[["qvalue", "pg_qvalue", "global_qvalue"]].apply(
            lambda row: [
                {"score_name": "qvalue", "score_value": row["qvalue"]},
                {"score_name": "pg_qvalue", "score_value": row["pg_qvalue"]},
                {"score_name": "global_qvalue", "score_value": row["global_qvalue"]},
            ],
            axis=1,
        )
        report.loc[:, "cv_params"] = report[["precursor_quantification_score"]].apply(
            lambda rows: [
                {"cv_name": "precursor_quantification_score", "cv_value": str(rows["precursor_quantification_score"])}
            ],
            axis=1,
        )
        report.loc[:, "scan_reference_file_name"] = None
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

    def write_pg_matrix_to_file(self, output_path: str, file_num=20):
        info_list = self.get_unique_references("Run")
        info_list = [info_list[i : i + file_num] for i in range(0, len(info_list), file_num)]
        pqwriter = None
        for refs in info_list:
            report = self.get_report_from_database(refs, DIANN_PG_SQL)
            report.rename(columns=DIANN_PG_MAP, inplace=True)
            report.dropna(subset=["pg_accessions"], inplace=True)
            for ref in refs:
                df = report[report["reference_file_name"] == ref].copy()
                df = self.generate_pg_matrix(df)
                pg_parquet = pa.Table.from_pandas(df, schema=PG_SCHEMA)
                if not pqwriter:
                    pqwriter = pq.ParquetWriter(output_path, pg_parquet.schema)
                pqwriter.write_table(pg_parquet)
        close_file(pqwriter=pqwriter)
        self.destroy_duckdb_database()

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
        file_num: int = 50,
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

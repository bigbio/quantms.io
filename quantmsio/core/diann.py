import duckdb
import time
import logging
import re
import numpy as np
import pandas as pd
import os
import pyarrow.parquet as pq
import concurrent.futures
from pathlib import Path
from pyopenms import AASequence
from pyopenms.Constants import PROTON_MASS_U
from quantmsio.core.project import create_uuid_filename
from quantmsio.core.sdrf import SDRFHandler
from quantmsio.core.psm import Psm
from quantmsio.core.feature import Feature
from quantmsio.utils.pride_utils import get_peptidoform_proforma_version_in_mztab,generate_scan_number
from quantmsio.core.common import DIANN_MAP, QUANTMSIO_VERSION

MODIFICATION_PATTERN = re.compile(r"\((.*?)\)")


def get_exp_design_dfs(exp_design_file):
    with open(exp_design_file, "r") as f:
        data = f.readlines()
        empty_row = data.index("\n")
        f_table = [i.replace("\n", "").split("\t") for i in data[1:empty_row]]
        f_header = data[0].replace("\n", "").split("\t")
        f_table = pd.DataFrame(f_table, columns=f_header)
        f_table.loc[:, "run"] = f_table.apply(lambda x: _true_stem(x["Spectra_Filepath"]), axis=1)
        f_table.rename(columns={"Fraction_Group": "ms_run", "run": "Run"}, inplace=True)
        s_table = [i.replace("\n", "").split("\t") for i in data[empty_row + 1 :]][1:]
        s_header = data[empty_row + 1].replace("\n", "").split("\t")
        s_data_frame = pd.DataFrame(s_table, columns=s_header)

    return s_data_frame, f_table


def _true_stem(x):
    """
    Return the true stem of a file name, i.e. the
    file name without the extension.

    :param x: The file name
    :type x: str
    :return: The true stem of the file name
    :rtype: str

    Examples:
    >>> _true_stem("foo.mzML")
    'foo'
    >>> _true_stem("foo.d.tar")
    'foo'

    These examples can be tested with pytest:
    $ pytest -v --doctest-modules
    """
    stem = os.path.splitext(os.path.basename(x))[0]

    return stem


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


class DiaNNConvert:

    def __init__(self, diann_report, sdrf_path, duckdb_max_memory="16GB", duckdb_threads=4):
        self._sdrf_path = sdrf_path
        self._diann_path = diann_report
        self._modifications = SDRFHandler(sdrf_path).get_mods_dict()
        self._duckdb_name = create_uuid_filename("report-duckdb", ".db")
        self._duckdb = self.create_duckdb_from_diann_report(duckdb_max_memory, duckdb_threads)

    def create_duckdb_from_diann_report(self, max_memory, worker_threads):
        """
        This function creates a duckdb database from a diann report for fast performance queries. The database
        is created from the tab delimited format of diann and can handle really large datasets.
        :param report_path: The path to the diann report
        :return: A duckdb database
        """
        s = time.time()

        database = duckdb.connect(self._duckdb_name)
        database.execute("SET enable_progress_bar=true")

        if max_memory is not None:
            database.execute("SET max_memory='{}'".format(max_memory))
        if worker_threads is not None:
            database.execute("SET worker_threads='{}'".format(worker_threads))

        msg = database.execute("SELECT * FROM duckdb_settings() where name in ('worker_threads', 'max_memory')").df()
        logging.info("duckdb uses {} threads.".format(str(msg["value"][0])))
        logging.info("duckdb uses {} of memory.".format(str(msg["value"][1])))

        database.execute("CREATE TABLE diann_report AS SELECT * FROM '{}'".format(self._diann_path))
        database.execute("""CREATE INDEX idx_precursor_q ON diann_report ("Precursor.Id", "Q.Value")""")
        database.execute("""CREATE INDEX idx_run ON diann_report ("Run")""")

        et = time.time() - s
        logging.info("Time to create duckdb database {} seconds".format(et))
        return database

    def get_report_from_database(self, runs: list) -> pd.DataFrame:
        """
        This function loads the report from the duckdb database for a group of ms_runs.
        :param runs: A list of ms_runs
        :return: The report
        """
        s = time.time()
        database = self._duckdb.query(
            """
            select "File.Name", "Run", "Protein.Ids", "Genes", "RT", "Predicted.RT", "RT.Start", "RT.Stop", "Precursor.Id", "Q.Value", "Global.Q.Value", "PEP",
                   "Global.PG.Q.Value", "Modified.Sequence", "Stripped.Sequence", "Precursor.Charge", "Precursor.Normalised"
                    from diann_report
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
            select DISTINCT "Modified.Sequence" from diann_report
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
            FROM diann_report
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

    def main_report_df(
        self,
        qvalue_threshold: float,
        mzml_info_folder: str,
        file_num: int,
    ):
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
            group.sort_values(by="RT.Start", inplace=True)
            target.rename(
                columns={
                    "Retention_Time": "RT.Start",
                    "SpectrumID": "opt_global_spectrum_reference",
                    "Exp_Mass_To_Charge": "exp_mass_to_charge",
                },
                inplace=True,
            )
            target["RT.Start"] = target["RT.Start"] / 60
            res = pd.merge_asof(group, target, on="RT.Start", direction="nearest")
            return res

        # query duckdb
        best_ref_map = self.get_peptide_map_from_database()
        masses_map, modifications_map = self.get_masses_and_modifications_map()

        info_list = [
            mzml.replace("_mzml_info.tsv", "")
            for mzml in os.listdir(mzml_info_folder)
            if mzml.endswith("_mzml_info.tsv")
        ]
        info_list = [info_list[i : i + file_num] for i in range(0, len(info_list), file_num)]
        for refs in info_list:
            report = self.get_report_from_database(refs)
            # restrict
            report = report[report["Q.Value"] < qvalue_threshold]
            usecols = report.columns.to_list() + [
                "opt_global_spectrum_reference",
                "exp_mass_to_charge",
            ]
            with concurrent.futures.ThreadPoolExecutor(100) as executor:
                results = executor.map(intergrate_msg, refs)
            report = np.vstack([result.values for result in results])
            report = pd.DataFrame(report, columns=usecols)

            # cal value and mod
            mass_vector = report["Modified.Sequence"].map(masses_map)
            report["Calculate.Precursor.Mz"] = (
                mass_vector + (PROTON_MASS_U * report["Precursor.Charge"].values)
            ) / report["Precursor.Charge"].values
            report["modifications"] = report["Modified.Sequence"].apply(lambda x: find_modification(x))
            report["Modified.Sequence"] = report["Modified.Sequence"].map(modifications_map)
            # pep
            report["psm_reference_file_name"] = report["Precursor.Id"].map(best_ref_map)
            report["psm_scan_number"] = None
            report.rename(columns=DIANN_MAP, inplace=True)
            # add extra msg
            report = self.add_additional_msg(report)
            yield report

    def add_additional_msg(self, report: pd.DataFrame) -> pd.DataFrame:
        """
        Perform some transformations in the report dataframe to help with the generation of the psm and feature files.
        :param report: The report dataframe
        :return: The report dataframe with the transformations
        """
        report["reference_file_name"] = report["reference_file_name"].apply(lambda x: x.split(".")[0])

        report.loc[:, "is_decoy"] = "0"
        report.loc[:, "unique"] = report["pg_accessions"].apply(lambda x: "0" if ";" in str(x) else "1")

        report.loc[:, "pg_positions"] = None

        report["peptidoform"] = report[["sequence", "modifications"]].apply(
            lambda row: get_peptidoform_proforma_version_in_mztab(
                row["sequence"], row["modifications"], self._modifications
            ),
            axis=1,
        )
        report["scan_number"] = report["scan_number"].apply(generate_scan_number)
        report.loc[:,'gg_names'] = report['gg_names'].str.split(',')
        report.loc[:,"additional_scores"] = None
        report.loc[:,"modification_details"] = None
        report.loc[:,"cv_params"] = None
        report.loc[:,"quantmsio_version"] = QUANTMSIO_VERSION
        report.loc[:,"gg_accessions"] = None
        return report

    def generate_psm_file(self, report, psm_pqwriter, psm_output_path):
        psm = report.copy()
        psm.loc[:, "intensity_array"] = None
        psm.loc[:, "ion_mobility"] = None
        psm.loc[:, "mz_array"] = None
        psm.loc[:, "num_peaks"] = None
        psm.loc[:, "rank"] = None
        parquet_table = Psm.convert_to_parquet(psm, self._modifications)
        if not psm_pqwriter:
            psm_pqwriter = pq.ParquetWriter(psm_output_path, parquet_table.schema)
        psm_pqwriter.write_table(parquet_table)
        return psm_pqwriter

    def generate_psm_and_feature_file(
        self,
        qvalue_threshold: float,
        mzml_info_folder: str,
        design_file: str,
        psm_output_path: str,
        feature_output_path: str,
        file_num: int = 50,
    ):
        psm_pqwriter = None
        feature_pqwriter = None

        s_data_frame, f_table = get_exp_design_dfs(design_file)
        for report in self.main_report_df(qvalue_threshold, mzml_info_folder, file_num):
            s = time.time()
            psm_pqwriter = self.generate_psm_file(report, psm_pqwriter, psm_output_path)
            feature_pqwriter = self.generate_feature_file(
                report,
                s_data_frame,
                f_table,
                feature_pqwriter,
                feature_output_path,
            )
            et = time.time() - s
            logging.info("Time to generate psm and feature file {} seconds".format(et))
        if psm_pqwriter:
            psm_pqwriter.close()
        if feature_pqwriter:
            feature_pqwriter.close()

        self.destroy_duckdb_database()

    def destroy_duckdb_database(self):
        """
        This function destroys the duckdb database, closing connection and removeing it from the filesystem.
        """
        if self._duckdb_name and self._duckdb:
            self._duckdb.close()
            os.remove(self._duckdb_name)
            self._duckdb_name = None
            self._duckdb = None

    def generate_feature_file(
        self,
        report,
        s_data_frame,
        f_table,
        feature_pqwriter,
        feature_output_path,
    ):
        sdrf = pd.read_csv(
            self._sdrf_path,
            sep="\t",
            usecols=[
                "source name",
                "comment[data file]",
                "comment[technical replicate]",
            ],
        )
        samples = sdrf["source name"].unique()
        mixed_map = dict(zip(samples, range(1, len(samples) + 1)))
        sdrf["comment[data file]"] = sdrf["comment[data file]"].apply(lambda x: x.split(".")[0])
        sdrf = sdrf.set_index(["comment[data file]"])
        sdrf_map = sdrf.to_dict()["source name"]
        tec_map = sdrf.to_dict()["comment[technical replicate]"]
        report = report[report["intensity"] != 0]
        report = report.merge(
            (
                s_data_frame[["Sample", "MSstats_Condition", "MSstats_BioReplicate"]]
                .merge(f_table[["Fraction", "Sample", "Run"]], on="Sample")
                .rename(
                    columns={
                        "MSstats_BioReplicate": "BioReplicate",
                        "MSstats_Condition": "Condition",
                    }
                )
                .drop(columns=["Sample"])
            ),
            on="Run",
            validate="many_to_one",
        )
        report.rename(
            columns={
                "Condition": "condition",
                "Fraction": "fraction",
                "BioReplicate": "biological_replicate",
                "Run": "run",
            },
            inplace=True,
        )
        report.loc[:, "channel"] = "LFQ"
        report.loc[:, "sample_accession"] = report["reference_file_name"].map(sdrf_map)
        report.loc[:, "comment[technical replicate]"] = report["reference_file_name"].map(tec_map)
        report.loc[:, "run"] = report[["sample_accession", "comment[technical replicate]", "fraction"]].apply(
            lambda row: str(mixed_map[row["sample_accession"]])
            + "_"
            + str(row["comment[technical replicate]"])
            + "_"
            + str(row["fraction"]),
            axis=1,
        )
        report.drop(["comment[technical replicate]"], axis=1, inplace=True)
        parquet_table = Feature.convert_to_parquet(report, self._modifications)

        if not feature_pqwriter:
            feature_pqwriter = pq.ParquetWriter(feature_output_path, parquet_table.schema)
        feature_pqwriter.write_table(parquet_table)

        return feature_pqwriter

import logging
import time

import pandas as pd
import os
from pyopenms import AASequence, ModificationsDB
from pyopenms.Constants import PROTON_MASS_U
from pathlib import Path
import numpy as np
import re
from quantms_io.core.mztab import fetch_modifications_from_mztab_line
from quantms_io.core.project import create_uuid_filename
from quantms_io.utils.pride_utils import (get_peptidoform_proforma_version_in_mztab)
import pyarrow as pa
import pyarrow.parquet as pq
from quantms_io.core.feature import FeatureHandler
from quantms_io.core.feature_in_memory import FeatureInMemory
from quantms_io.core.psm import PSMHandler

import concurrent.futures
import duckdb
import swifter
MODIFICATION_PATTERN = re.compile(r"\((.*?)\)")


def get_exp_design_dfs(exp_design_file):
    # logger.info(f"Reading experimental design file: {exp_design_file}")
    with open(exp_design_file, "r") as f:
        data = f.readlines()
        empty_row = data.index("\n")
        f_table = [i.replace("\n", "").split("\t") for i in data[1:empty_row]]
        f_header = data[0].replace("\n", "").split("\t")
        f_table = pd.DataFrame(f_table, columns=f_header)
        f_table.loc[:, "run"] = f_table.apply(lambda x: _true_stem(x["Spectra_Filepath"]), axis=1)
        f_table.rename(columns={"Fraction_Group": "ms_run", "run": "Run"}, inplace=True)

        s_table = [i.replace("\n", "").split("\t") for i in data[empty_row + 1:]][1:]
        s_header = data[empty_row + 1].replace("\n", "").split("\t")
        s_DataFrame = pd.DataFrame(s_table, columns=s_header)

    return s_DataFrame, f_table


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

def mtd_mod_info(fix_mod, var_mod):
    """
    Convert fixed and variable modifications to the format required by the MTD sub-table.

    :param fix_mod: Fixed modifications from DIA parameter list
    :type fix_mod: str
    :param var_mod: Variable modifications from DIA parameter list
    :type var_mod: str
    :return: A tuple contains fixed and variable modifications, and flags indicating whether they are null
    :rtype: tuple
    """
    var_ptm = []
    fix_ptm = []
    mods_db = ModificationsDB()

    if fix_mod != "null":
        fix_flag = 1
        for mod in fix_mod.split(","):
            mod_obj = mods_db.getModification(mod)
            mod_name = mod_obj.getId()
            mod_accession = mod_obj.getUniModAccession()
            site = mod_obj.getOrigin()
            fix_ptm.append(("[UNIMOD, " + mod_accession.upper() + ", " + mod_name + ", ]", site))
    else:
        fix_flag = 0
        fix_ptm.append("[MS, MS:1002453, No fixed modifications searched, ]")

    if var_mod != "null":
        var_flag = 1
        for mod in var_mod.split(","):
            mod_obj = mods_db.getModification(mod)
            mod_name = mod_obj.getId()
            mod_accession = mod_obj.getUniModAccession()
            site = mod_obj.getOrigin()
            var_ptm.append(("[UNIMOD, " + mod_accession.upper() + ", " + mod_name + ", ]", site))
    else:
        var_flag = 0
        var_ptm.append("[MS, MS:1002454, No variable modifications searched, ]")

    return fix_ptm, var_ptm, fix_flag, var_flag


def get_modifications(fix_modifications: str, variable_modifications: str):
    (fixed_mods, variable_mods, fix_flag, var_flag) = mtd_mod_info(fix_modifications, variable_modifications)

    modifications_list = []
    if fix_flag == 1:
        for i in range(1, len(fixed_mods) + 1):
            modifications_list.append(f"MTD\tfixed_mod[{str(i)}]\t{fixed_mods[i - 1][0]}")
            modifications_list.append(f"MTD\tfixed_mod[{str(i)}]-site\t{fixed_mods[i - 1][1]}")
            modifications_list.append(f"MTD\tfixed_mod[{str(i)}]-position\tAnywhere")
    else:
        modifications_list.append(f"MTD\tfixed_mod[1]\t{fixed_mods[0]}")

    if var_flag == 1:
        for i in range(1, len(variable_mods) + 1):
            modifications_list.append(f"MTD\tvariable_mod[{str(i)}]\t{variable_mods[i - 1][0]}")
            modifications_list.append(f"MTD\tvariable_mod[{str(i)}]-site\t{variable_mods[i - 1][1]}")
            modifications_list.append(f"MTD\tvariable_mod[{str(i)}]-position\tAnywhere")
    else:
        modifications_list.append(f"MTD\tvariable_mod[1]\t{variable_mods[0]}")

    mod_dict = {}
    for line in modifications_list:
        mod_dict = fetch_modifications_from_mztab_line(line, mod_dict)

    return mod_dict



class DiaNNConvert:

    def __init__(self):

        self._modifications = None
        self._score_names = {'peptide_score': 'DIA-NN Q.Value (minimum of the respective precursor q-values)',
                             'protein_score': 'DIA-NN Global.PG.Q.Value',
                             'psm_score': 'protein-level q-value'}
        self._columns_map = {
                            "Stripped.Sequence": "sequence",
                            "Protein.Ids": "protein_accessions",
                            "RT.Start": "retention_time",
                            "Precursor.Charge": "charge",
                            "Calculate.Precursor.Mz": "calc_mass_to_charge",
                            "exp_mass_to_charge": "exp_mass_to_charge",
                            "Modified.Sequence": "peptidoform",
                            "opt_global_spectrum_reference": "scan_number",
                            "File.Name": "reference_file_name",
                            "Precursor.Quantity": "intensity",
                            "Global.Q.Value": "global_qvalue",
                            "PEP": "posterior_error_probability",
                            "Global.PG.Q.Value": "protein_global_qvalue",
        }

    def create_duckdb_from_diann_report(self, report_path, max_memory, worker_threads):
        """
        This function creates a duckdb database from a diann report for fast performance queries. The database
        is created from the tab delimited format of diann and can handle really large datasets.
        :param report_path: The path to the diann report
        :return: A duckdb database
        """
        s = time.time()

        self._duckdb_name = create_uuid_filename("report-duckdb", '.db')

        database = duckdb.connect(self._duckdb_name)
        database.execute("SET enable_progress_bar=true")

        if max_memory is not None:
            database.execute("SET max_memory='{}'".format(max_memory))
        if worker_threads is not None:
            database.execute("SET worker_threads='{}'".format(worker_threads))

        msg = database.execute("SELECT * FROM duckdb_settings() where name in ('worker_threads', 'max_memory')").df()
        logging.info('Duckdb uses {} threads.'.format(str(msg['value'][0])))
        logging.info('duckdb uses {} of memory.'.format(str(msg['value'][1])))

        database.execute("CREATE TABLE diann_report AS SELECT * FROM '{}'".format(report_path))
        database.execute("""CREATE INDEX idx_precursor_q ON diann_report ("Precursor.Id", "Q.Value")""")
        database.execute("""CREATE INDEX idx_run ON diann_report ("Run")""")

        et = time.time() - s
        logging.info('Time to create duckdb database {} seconds'.format(et))
        return database

    def get_protein_map_from_database(self):
        s = time.time()
        query = self._duckdb.query(
            """
            select "Protein.Ids",MIN("Global.PG.Q.Value") as "Global.PG.Q.Value" from diann_report
            GROUP BY "Protein.Ids"
            """
        )
        protein_df = query.df()
        protein_df.index = protein_df["Protein.Ids"]
        protein_map = protein_df.to_dict()["Global.PG.Q.Value"]

        et = time.time() - s
        logging.info('Time to load protein map {} seconds'.format(et))

        return protein_map

    def get_peptide_map_from_database(self):
        s = time.time()
        database = self._duckdb.query(
            """    
            SELECT "Precursor.Id","Q.Value","Run"
            FROM (
            SELECT "Precursor.Id", "Q.Value","Run", ROW_NUMBER() OVER (PARTITION BY "Precursor.Id" ORDER BY "Q.Value" ASC) AS row_num
            FROM diann_report
            ) AS subquery
            WHERE row_num = 1;
            """
        )
        peptide_df = database.df()
        peptide_df.set_index("Precursor.Id", inplace=True)
        peptide_map = peptide_df.to_dict()["Q.Value"]
        best_ref_map = peptide_df.to_dict()["Run"]
        et = time.time()-s
        logging.info('Time to load peptide map {} seconds'.format(et))
        return peptide_map, best_ref_map
    
    def get_report_from_database(self, runs: list) -> pd.DataFrame:
        """
        This function loads the report from the duckdb database for a group of ms_runs.
        :param runs: A list of ms_runs
        :return: The report
        """
        s = time.time()
        database = self._duckdb.query(
            """
            select "File.Name", "Run", "Protein.Ids", "RT.Start", "Precursor.Id", "Q.Value", "Global.Q.Value", "PEP", 
                   "Global.PG.Q.Value", "Modified.Sequence", "Stripped.Sequence", "Precursor.Charge", "Precursor.Quantity"
                    from diann_report
            where Run IN {}
            """.format(tuple(runs))
        )
        report = database.df()
        et = time.time() - s
        logging.info('Time to load report {} seconds'.format(et))
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
        
        return masses_map,modifications_map
    
    def main_report_df(self, report_path: str, qvalue_threshold: float, mzml_info_folder: str, file_num: int) -> pd.DataFrame:
        def intergrate_msg(n):
            nonlocal report
            nonlocal mzml_info_folder
            files = list(Path(mzml_info_folder).glob(f"*{n}_mzml_info.tsv"))
            if not files:
                raise ValueError(f"Could not find {n} info file in {dir}")
            target = pd.read_csv(files[0],sep='\t',usecols=["Retention_Time", "SpectrumID", "Exp_Mass_To_Charge"])
            group = report[report['Run']==n].copy()
            group.sort_values(by="RT.Start", inplace=True)
            target.rename(columns={"Retention_Time": "RT.Start", "SpectrumID": "opt_global_spectrum_reference",
                                        "Exp_Mass_To_Charge": "exp_mass_to_charge"}, inplace=True)
            target["RT.Start"] = target["RT.Start"] / 60
            res = pd.merge_asof(group, target, on="RT.Start", direction="nearest")
            return res

        #query duckdb
        peptide_map, best_ref_map = self.get_peptide_map_from_database()
        masses_map, modifications_map = self.get_masses_and_modifications_map()
        
        info_list = [mzml.replace('_mzml_info.tsv','') for mzml in os.listdir(mzml_info_folder) if mzml.endswith('_mzml_info.tsv')]
        info_list =  [info_list[i:i+file_num] for i in range(0,len(info_list), file_num)]
        for refs in info_list:
            report = self.get_report_from_database(refs)

            #restrict
            report = report[report["Q.Value"] < qvalue_threshold]
            usecols = report.columns.to_list() + ["opt_global_spectrum_reference", "exp_mass_to_charge"]
            with concurrent.futures.ThreadPoolExecutor(100) as executor:
                results = executor.map(intergrate_msg,refs)
            report = np.vstack([result.values for result in results])
            report = pd.DataFrame(report, columns=usecols)

            #spectral count
            report["spectral_count"] = report.groupby(["Modified.Sequence","Precursor.Charge","Run"]).transform('size')
            #cal value and mod
            mass_vector = report["Modified.Sequence"].map(masses_map)
            report["Calculate.Precursor.Mz"] = (mass_vector + (PROTON_MASS_U * report["Precursor.Charge"].values)) / report["Precursor.Charge"].values
            report["modifications"] = report["Modified.Sequence"].swifter.apply(lambda x: find_modification(x))
            report["Modified.Sequence"] = report["Modified.Sequence"].map(modifications_map)
            #pep 
            report["best_psm_reference_file_name"] = report["Precursor.Id"].map(best_ref_map)
            report["best_psm_scan_number"] = None
            # add extra msg
            report = self.add_additional_msg(report)
            yield report

    def generate_psm_and_feature_file(self, report_path: str, qvalue_threshold: float, mzml_info_folder: str,
                                      design_file: str, modifications:list, sdrf_path:str, psm_output_path:str,
                                      feature_output_path:str, duckdb_max_memory:str = None, duckdb_threads:int = None, file_num:int=2):
        psm_pqwriter = None
        feature_pqwriter = None

        self._duckdb = self.create_duckdb_from_diann_report(report_path, duckdb_max_memory, duckdb_threads)

        s_data_frame, f_table = get_exp_design_dfs(design_file)
        self._modifications = get_modifications(modifications[0], modifications[1])
        for report in self.main_report_df(report_path, qvalue_threshold, mzml_info_folder , file_num):
            s = time.time()
            psm_pqwriter = self.generate_psm_file(report, psm_pqwriter, psm_output_path)
            feature_pqwriter = self.generate_feature_file(report,s_data_frame,f_table,sdrf_path,feature_pqwriter,feature_output_path)
            et = time.time() - s
            logging.info('Time to generate psm and feature file {} seconds'.format(et))
        if psm_pqwriter:
            psm_pqwriter.close()
        if feature_pqwriter:
            feature_pqwriter.close()

        self.destroy_duckdb_database()

    def add_additional_msg(self,report:pd.DataFrame) -> pd.DataFrame:
        """
        Perform some transformations in the report dataframe to help with the generation of the psm and feature files.
        :param report: The report dataframe
        :return: The report dataframe with the transformations
        """
        report.rename(columns=self._columns_map , inplace=True)
        report["reference_file_name"] = report["reference_file_name"].swifter.apply(lambda x:x.split('.')[0])

        report.loc[:,"is_decoy"] = "0"
        report.loc[:, "unique"] = report["protein_accessions"].swifter.apply(lambda x: "0" if ";" in str(x) else "1")

        null_col = [
            "protein_start_positions",
            "protein_end_positions"
        ]
        report.loc[:, null_col] = None

        report['peptidoform'] = report[['sequence', 'modifications']].swifter.apply(lambda row:
            get_peptidoform_proforma_version_in_mztab(row['sequence'], row['modifications'], self._modifications),
                                                                                    axis=1)

        report['id_scores'] = report[['Q.Value', 'posterior_error_probability', 'global_qvalue']].swifter.apply(lambda x: [
            f"q-value: {x['Q.Value']}",
            f"global q-value: {x['global_qvalue']}",
            f"posterior error probability: {x['posterior_error_probability']}",
        ], axis=1)
        
        return report
    
    def generate_psm_file(self,report,psm_pqwriter,psm_output_path):
        psm = report.copy()
        psm.loc[:, "consensus_support"] = None
        parquet_table = self.convert_psm_format_to_parquet(psm)
        
        if not psm_pqwriter:
            # create a parquet write object giving it an output file
            psm_pqwriter = pq.ParquetWriter(psm_output_path, parquet_table.schema)
        psm_pqwriter.write_table(parquet_table)
        return psm_pqwriter

    def generate_feature_file(self,report,s_data_frame,f_table,sdrf_path,feature_pqwriter,feature_output_path):

        sdrf = pd.read_csv(sdrf_path,sep='\t',usecols=['source name',"comment[data file]"])
        sdrf["comment[data file]"] = sdrf["comment[data file]"].apply(lambda x: x.split('.')[0])
        sdrf = sdrf.set_index(["comment[data file]"])
        sdrf_map  = sdrf.to_dict()['source name']
        report = report[report["intensity"] != 0]
        report.loc[:,"fragment_ion"] = "NA"
        report.loc[:,"isotope_label_type"] = "L"
        report.loc[:,"channel"] = "LABEL FREE SAMPLE"
        report = report.merge(
            (
                s_data_frame[["Sample", "MSstats_Condition", "MSstats_BioReplicate"]]
                .merge(f_table[["Fraction", "Sample", "Run"]], on="Sample")
                .rename(
                    columns={"MSstats_BioReplicate": "BioReplicate", "MSstats_Condition": "Condition"})
                .drop(columns=["Sample"])
            ),
            on="Run",
            validate="many_to_one",
        )
        report.rename(columns={
            'Condition': 'condition',
            'Fraction': 'fraction',
            'BioReplicate': 'biological_replicate',
            'Run': 'run'
        },inplace=True)
        peptide_score_name = self._score_names["peptide_score"]
        report['sample_accession'] = report["reference_file_name"].map(sdrf_map)
        schema = FeatureHandler()
        feature = FeatureInMemory('LFQ',schema.schema)
        feature._modifications = self._modifications
        feature._score_names = self._score_names
        parquet_table = feature.convert_to_parquet(report)

        if not feature_pqwriter:
            # create a parquet write object giving it an output file
            feature_pqwriter = pq.ParquetWriter(feature_output_path, parquet_table.schema)
        feature_pqwriter.write_table(parquet_table)
        
        return feature_pqwriter
    

    @staticmethod
    def __split_start_or_end(value):
        """
        split start or end
        :param value: start or end
        """
        if pd.isna(value) or value == 'null':
            return pd.NA
        elif "," in str(value):
            return list(map(int, value.split(",")))
        elif value is np.nan:
            return None
        else:
            return [int(value)]
    
    def convert_psm_format_to_parquet(self, res):
        """
            res: msstats_in dataframe
            return: parquet table
        """
        psm = PSMHandler()
        feature = FeatureInMemory('LFQ', None)
        feature._modifications = self._modifications
        feature._score_names = self._score_names
        res['sequence'] = res['sequence'].astype(str)
        res['protein_accessions'] = res['protein_accessions'].str.split(";")
        res['protein_start_positions'] = res['protein_start_positions'].swifter.apply(
            self.__split_start_or_end).to_list()
        res['protein_end_positions'] = res['protein_end_positions'].swifter.apply(
            self.__split_start_or_end).to_list()
        res['protein_global_qvalue'] = res['protein_global_qvalue'].astype(float)
        res['unique'] = res['unique'].map(lambda x: pd.NA if pd.isna(x) else int(x)).astype('Int32')
        res['modifications'] = res['modifications'].swifter.apply(lambda x: feature._generate_modification_list(x))
        res['retention_time'] = res['retention_time'].astype(float)
        res['charge'] = res['charge'].map(lambda x: pd.NA if pd.isna(x) else int(x)).astype('Int32')
        res['exp_mass_to_charge'] = res['exp_mass_to_charge'].astype(float)
        res['calc_mass_to_charge'] = res['calc_mass_to_charge'].astype(float)
        res['posterior_error_probability'] = res['posterior_error_probability'].astype(float)
        res['global_qvalue'] = res['global_qvalue'].astype(float)

        res['is_decoy'] = res['is_decoy'].map(lambda x: pd.NA if pd.isna(x) else int(x)).astype('Int32')

        res.loc[:, "num_peaks"] = None
        res.loc[:, "mz_array"] = None
        res.loc[:, "intensity_array"] = None

        res.loc[:, "gene_accessions"] = None
        res.loc[:, "gene_names"] = None

        return pa.Table.from_pandas(res, schema=psm._create_schema())

    def destroy_duckdb_database(self):
        """
        This function destroys the duckdb database, closing connection and removeing it from the filesystem.
        """
        if self._duckdb_name and self._duckdb:
            self._duckdb.close()
            os.remove(self._duckdb_name)
            self._duckdb_name = None
            self._duckdb = None

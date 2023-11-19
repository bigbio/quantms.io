import pandas as pd
import os
import time
from pyopenms import AASequence, FASTAFile, ModificationsDB
from pyopenms.Constants import PROTON_MASS_U
from pathlib import Path
import numpy as np
import re
from quantms_io.core.mztab import fetch_modifications_from_mztab_line
from quantms_io.utils.pride_utils import (clean_peptidoform_sequence, get_peptidoform_proforma_version_in_mztab,
                                          generate_scan_number, print_estimated_time)
import pyarrow as pa
import pyarrow.parquet as pq
from quantms_io.core.feature import FeatureHandler
from quantms_io.core.feature_in_memory import FeatureInMemory
from quantms_io.core.psm import PSMHandler
from collections import Counter
from quantms_io.utils.thread import MyThread
from typing import Any, List, Tuple, Dict, Set

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


def generate_scan_number(spectra_ref: str):
    if 'scan' in spectra_ref:
        return re.findall(r"scan=(\d+)", spectra_ref)[0]
    else:
        return ",".join(re.findall(r'=(\d+)', spectra_ref))


def handle_protein_map(protein_map, key):
    """
    map protein score from accession
    """
    key = key.replace(";", ",")
    if key not in protein_map.keys():
        keys = key.split(";")
        for k in keys:
            if k in protein_map.keys():
                return protein_map[k]
        return None
    else:
        return protein_map[key]


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

        self._ms_runs = None
        self._modifications = None
        self._protein_map = None
        self.uniq_masses_map = None
        self.unique_reference_map = None
        self.protein_best_score_dict = None
        self.peptide_best_score_dict = None
        self._score_names = {'peptide_score': 'DIA-NN Q.Value (minimum of the respective precursor q-values)',
                             'protein_score': 'DIA-NN Global.PG.Q.Value',
                             'psm_score': 'protein-level q-value'}
        self._lfq_msstats_usecols = [
            "ProteinName",
            "PeptideSequence",
            "Channel",
            "Reference",
            "Run",
            "BioReplicate",
            "Intensity",
            "FragmentIon",
            "IsotopeLabelType",
            "PrecursorCharge",
            "protein_global_qvalue",
            "sequence",
        ]
        self.map_psm = {
            "accession": "protein_accessions",
            "start": "protein_start_positions",
            "end": "protein_end_positions",
            "spectra_ref": "reference_file_name",
            "opt_global_q-value": "global_qvalue",
            "opt_global_cv_MS:1002217_decoy_peptide": "is_decoy"
        }

    def get_uniq_and_index_and_reference_map(self, report_path, chunksize: int):
        use_cols = ["Modified.Sequence", "Precursor.Id", "File.Name", "Protein.Ids", "Global.PG.Q.Value", "Q.Value"]
        reports = pd.read_csv(report_path, sep="\t", header=0, usecols=use_cols, chunksize=chunksize)
        uniq_set = set()
        # precursor_indexs = set()
        reference_set = set()

        protein_best_score_dict = {}
        peptide_best_score_dict = {}

        for report in reports:
            uniq_set.update(report["Modified.Sequence"].unique())
            # precursor_indexs.update(report["Precursor.Id"].unique())
            reference_set.update(report["File.Name"].unique())

            grouped_df = (
                report[["Modified.Sequence", "Protein.Ids", "Global.PG.Q.Value"]]
                .sort_values("Global.PG.Q.Value", ascending=True)
                .groupby(["Protein.Ids"])
                .head(1)
            )
            out = {
                row["Protein.Ids"]: (row["Modified.Sequence"], row["Global.PG.Q.Value"]) for _, row in
                grouped_df.iterrows()
            }
            protein_best_score_dict = self.__update_dict(protein_best_score_dict, out, 1)

            aggtable = (report.groupby(["Precursor.Id"])
            .agg(
                {
                    "Q.Value": "min",
                }
            )
            .reset_index()
            .rename(
                columns={
                    "Q.Value": "best_search_engine_score[1]",
                }
            )
            )
            aggtable.index = aggtable["Precursor.Id"]
            pr_dict = aggtable.to_dict()['best_search_engine_score[1]']
            peptide_best_score_dict = self.__update_dict(peptide_best_score_dict, pr_dict, 0)

        uniq_masses_map = {k: AASequence.fromString(k).getMonoWeight() for k in uniq_set}
        # precursor_index_map = {k: i for i, k in enumerate(list(precursor_indexs))}
        unique_reference_map = {k: os.path.basename(k) for k in reference_set}

        return uniq_masses_map, unique_reference_map, protein_best_score_dict, peptide_best_score_dict

    def __update_dict(self, map_dict, temporary_dict, n):
        if len(map_dict) == 0:
            map_dict.update(temporary_dict)
        else:
            for key in temporary_dict.keys():
                if key in map_dict:
                    if n != 0:
                        if float(map_dict[key][n]) > float(temporary_dict[key][n]):
                            map_dict[key] = temporary_dict[key]
                    else:
                        if float(map_dict[key]) > float(temporary_dict[key]):
                            map_dict[key] = temporary_dict[key]
                else:
                    map_dict[key] = temporary_dict[key]
        return map_dict

    def main_report_df(self, report_path: str, qvalue_threshold: float, chunksize: int,
                       uniq_masses: dict) -> pd.DataFrame:
        """
        The main_report_df method is a part of the DiaNNConvert class and is used to process a report file and filter
        the data based on a given q-value threshold. It returns a pandas DataFrame containing the filtered data.

        :param report_path (str): The path to the report file.
        :param qvalue_threshold (float): The q-value threshold for filtering the data.
        :param chunksize (int): The number of rows to read from the report file at a time.
        :param uniq_masses (dict): A dictionary containing unique masses for modified sequences.
        :return (pd.DataFrame): A pandas DataFrame containing the filtered data.
        """

        remain_cols = [
            "File.Name",
            "Run",
            "Protein.Group",
            "Protein.Names",
            "Protein.Ids",
            "First.Protein.Description",
            "PG.MaxLFQ",
            "RT.Start",
            "Global.Q.Value",
            "Lib.Q.Value",
            "PEP",
            "Precursor.Normalised",
            "Precursor.Id",
            "Q.Value",
            "Modified.Sequence",
            "Stripped.Sequence",
            "Precursor.Charge",
            "Precursor.Quantity",
            "Global.PG.Q.Value",
        ]
        reports = pd.read_csv(report_path, sep="\t", header=0, usecols=remain_cols, chunksize=chunksize)

        for report in reports:
            st = time.time()
            report = report[report["Q.Value"] < qvalue_threshold]
            mass_vector = report["Modified.Sequence"].map(uniq_masses)
            report["Calculate.Precursor.Mz"] = (mass_vector + (PROTON_MASS_U * report["Precursor.Charge"])) / report[
                "Precursor.Charge"
            ]
            # report["precursor.Index"] = report["Precursor.Id"].map(precursor_index_map)
            print_estimated_time(st, "{} create report chunksize".format(chunksize))
            yield report

    def get_msstats_in(self, report, unique_reference_map, s_data_frame, f_table):
        msstats_columns_keep = [
            "Protein.Names",
            "Modified.Sequence",
            "Precursor.Charge",
            "Precursor.Quantity",
            "File.Name",
            "Run",
        ]
        out_msstats = report[msstats_columns_keep]
        out_msstats.columns = ["ProteinName", "PeptideSequence", "PrecursorCharge", "Intensity", "Reference", "Run"]
        out_msstats = out_msstats[out_msstats["Intensity"] != 0]
        out_msstats.loc[:, "PeptideSequence"] = out_msstats.apply(
            lambda x: AASequence.fromString(x["PeptideSequence"]).toString(), axis=1)
        out_msstats["FragmentIon"] = "NA"
        out_msstats["ProductCharge"] = "0"
        out_msstats["IsotopeLabelType"] = "L"
        out_msstats["Reference"] = out_msstats["Reference"].map(unique_reference_map)
        out_msstats = out_msstats.merge(
            (
                s_data_frame[["Sample", "MSstats_Condition", "MSstats_BioReplicate"]]
                .merge(f_table[["Fraction", "Sample", "run"]], on="Sample")
                .rename(
                    columns={"run": "Run", "MSstats_BioReplicate": "BioReplicate", "MSstats_Condition": "Condition"})
                .drop(columns=["Sample"])
            ),
            on="Run",
            validate="many_to_one",
        )
        return out_msstats

    def generete_fasta_msg(self, fasta_path: str, f_table: pd.DataFrame):
        entries = []
        f = FASTAFile()
        f.load(fasta_path, entries)
        fasta_entries = [(e.identifier, e.sequence, len(e.sequence)) for e in entries]
        fasta_df = pd.DataFrame(fasta_entries, columns=["id", "seq", "len"])

        index_ref = f_table.copy()
        index_ref.rename(columns={"Fraction_Group": "ms_run", "Sample": "study_variable", "run": "Run"}, inplace=True)
        index_ref["ms_run"] = index_ref["ms_run"].astype("int")
        index_ref["study_variable"] = index_ref["study_variable"].astype("int")

        return fasta_df, index_ref

    def __map_ms_runs(self, ms_run, spectra, ms_runs):
        run = f"ms_run[{ms_run}]"
        ms_runs[run] = spectra.split('.')[0]

    def get_ms_runs(self, index_ref):
        ms_runs = {}
        index_ref[['ms_run', "Spectra_Filepath"]].apply(
            lambda row: self.__map_ms_runs(row['ms_run'], row["Spectra_Filepath"], ms_runs), axis=1)
        return ms_runs

    def get_protein_map(self, pg_path: str, best_score_dict: dict):
        pg = pd.read_csv(
            pg_path,
            sep="\t",
            header=0,
        )

        pg.loc[pg["Protein.Ids"].str.contains(";"), "opt_global_result_type"] = "indistinguishable_protein_group"

        out_mztab_prh = pg.drop(["Protein.Names"], axis=1)
        out_mztab_prh.rename(
            columns={"Protein.Group": "accession", "First.Protein.Description": "description"}, inplace=True
        )
        out_mztab_prh.loc[:, "accession"] = out_mztab_prh.apply(lambda x: x["accession"].split(";")[0], axis=1)

        protein_details_df = out_mztab_prh[out_mztab_prh["opt_global_result_type"] == "indistinguishable_protein_group"]
        prh_series = protein_details_df["Protein.Ids"].str.split(";", expand=True).stack().reset_index(level=1,
                                                                                                       drop=True)
        prh_series.name = "accession"
        protein_details_df = (
            protein_details_df.drop("accession", axis=1).join(prh_series).reset_index().drop(columns="index")
        )
        out_mztab_prh = pd.concat([out_mztab_prh, protein_details_df]).reset_index(drop=True)
        out_mztab_prh.loc[:, "ambiguity_members"] = out_mztab_prh.apply(
            lambda x: x["Protein.Ids"] if x["opt_global_result_type"] == "indistinguishable_protein_group" else "null",
            axis=1,
        )
        out_mztab_prh[["modifiedSequence", "best_search_engine_score[1]"]] = out_mztab_prh.apply(
            lambda x: best_score_dict.get(x["Protein.Ids"], (np.nan, np.nan)), axis=1, result_type="expand"
        )
        prt_score = out_mztab_prh[["ambiguity_members", "best_search_engine_score[1]"]].groupby(
            "ambiguity_members").min()
        protein_map = prt_score.to_dict()["best_search_engine_score[1]"]

        return protein_map

    def extract_dict_from_pep(self, pr_path: str, bset_score_map: dict) -> dict:
        pep_usecols = ['Protein.Group',
                       'Protein.Ids',
                       'Protein.Names',
                       'Genes',
                       'First.Protein.Description',
                       'Proteotypic',
                       'Stripped.Sequence',
                       'Modified.Sequence',
                       'Precursor.Charge',
                       'Precursor.Id']
        out_matab_peh = pd.read_csv(
            pr_path,
            sep="\t",
            header=0,
            usecols=pep_usecols)

        out_matab_peh.drop(
            ["Protein.Group", "Protein.Names", "First.Protein.Description", "Proteotypic"], axis=1, inplace=True
        )
        out_matab_peh.rename(
            columns={
                "Stripped.Sequence": "sequence",
                "Protein.Ids": "accession",
                "Modified.Sequence": "opt_global_cv_MS:1000889_peptidoform_sequence",
                "Precursor.Charge": "charge",
            },
            inplace=True,
        )
        out_matab_peh.loc[:, "opt_global_cv_MS:1000889_peptidoform_sequence"] = out_matab_peh.apply(
            lambda x: AASequence.fromString(x["opt_global_cv_MS:1000889_peptidoform_sequence"]).toString(), axis=1
        )

        out_matab_peh["best_search_engine_score[1]"] = out_matab_peh["Precursor.Id"].map(bset_score_map)
        out_matab_peh.fillna("null", inplace=True)
        out_matab_peh.loc[:, "scan_number"] = None
        out_matab_peh.loc[:, "spectra_ref"] = None
        pep = out_matab_peh[
            ['charge', 'opt_global_cv_MS:1000889_peptidoform_sequence', 'best_search_engine_score[1]', "scan_number",
             "spectra_ref"]]

        pep_msg = pep.iloc[
            pep.groupby(
                ["opt_global_cv_MS:1000889_peptidoform_sequence", "charge"]
            ).apply(lambda row: row["best_search_engine_score[1]"].idxmin())
        ]
        pep_msg = pep_msg.set_index(
            ["opt_global_cv_MS:1000889_peptidoform_sequence", "charge"]
        )

        pep_msg.loc[:, "pep_msg"] = pep_msg[
            ["best_search_engine_score[1]", "spectra_ref", "scan_number"]
        ].apply(
            lambda row: [
                row["best_search_engine_score[1]"],
                row["spectra_ref"],
                row["scan_number"],
            ],
            axis=1,
        )

        map_dict = pep_msg.to_dict()["pep_msg"]

        return map_dict

    def extract_from_psm_to_pep_msg(self, report_path: str, qvalue_threshold: float, folder: str,
                                    uniq_masses: Dict, index_ref: Dict, map_dict: Dict, psm_output_path:str, chunksize: int):
        """
        The extract_from_psm_to_pep_msg method is part of the DiaNNConvert class and is responsible for extracting
        relevant information from a PSM (Peptide-Spectrum Match) report and converting it into a peptide message format.

        :param report_path: The path to the PSM report file in TSV format.
        :param qvalue_threshold: The threshold value for the Q.Value score.
        :param folder: The path to the folder containing additional information files.
        :param uniq_masses: A dictionary mapping unique peptide sequences to their corresponding mass values.
        :param index_ref: A dictionary containing index and reference mappings.
        :param map_dict: A dictionary used for mapping values.
        :param chunksize: The number of rows to read at a time from the PSM report file.
        :return: A pandas DataFrame containing the extracted information.
        """

        psm_unique_keys = []
        spectra_count_dict = Counter()

        pqwriter = None
        for report in self.main_report_df(report_path, qvalue_threshold, chunksize, uniq_masses):

            report = report.merge(index_ref[["ms_run", "Run", "study_variable"]], on="Run", validate="many_to_one")


            def intergrate_msg(folder,n,group):
                files = list(Path(folder).glob(f"*{n}_mzml_info.tsv"))
                if not files:
                    raise ValueError(f"Could not find {n} info file in {dir}")
                target = pd.read_csv(files[0],sep='\t',usecols=["Retention_Time", "SpectrumID", "Exp_Mass_To_Charge"])
                group.sort_values(by="RT.Start", inplace=True)
                target.rename(columns={"Retention_Time": "RT.Start", "SpectrumID": "opt_global_spectrum_reference",
                                            "Exp_Mass_To_Charge": "exp_mass_to_charge"}, inplace=True)
                target["RT.Start"] = target["RT.Start"] / 60
                res = pd.merge_asof(group, target, on="RT.Start", direction="nearest")
                return res
            
            out_mztab_psh = pd.DataFrame()
            ThreadPool = []
            references = report['Run'].unique().tolist()
            references_threadpools =  [references[i:i+16] for i in range(0,len(references),16)]

            for refs in references_threadpools:
                for ref in refs:
                    ThreadPool.append(MyThread(target=intergrate_msg, args=(folder,ref,report[report['Run']==ref].copy())))
                for p in ThreadPool:
                    p.start()
                for p in ThreadPool:
                    p.join()
                out_mztab_psh = pd.concat([out_mztab_psh]+[p.result for p in ThreadPool])
                ThreadPool = []



            ## Score at PSM level: Q.Value
            out_mztab_psh = out_mztab_psh[
                [
                    "Stripped.Sequence",
                    "Protein.Ids",
                    "Q.Value",
                    "RT.Start",
                    "Precursor.Charge",
                    "Calculate.Precursor.Mz",
                    "exp_mass_to_charge",
                    "Modified.Sequence",
                    "PEP",
                    "Global.Q.Value",
                    "Global.Q.Value",
                    "opt_global_spectrum_reference",
                    "ms_run",
                ]
            ]
            out_mztab_psh.columns = [
                "sequence",
                "accession",
                "search_engine_score[1]",
                "retention_time",
                "charge",
                "calc_mass_to_charge",
                "exp_mass_to_charge",
                "opt_global_cv_MS:1000889_peptidoform_sequence",
                "opt_global_SpecEValue_score",
                "opt_global_q-value",
                "opt_global_q-value_score",
                "opt_global_spectrum_reference",
                "ms_run",
            ]

            out_mztab_psh.loc[:, "opt_global_cv_MS:1002217_decoy_peptide"] = "0"
            out_mztab_psh.loc[:, "PSM_ID"] = out_mztab_psh.index
            out_mztab_psh.loc[:, "unique"] = out_mztab_psh.apply(lambda x: "0" if ";" in str(x["accession"]) else "1",
                                                                 axis=1, result_type="expand")

            null_col = [
                "start",
                "end",
                "opt_global_feature_id",
                "opt_global_map_index",
            ]
            out_mztab_psh.loc[:, null_col] = "null"

            out_mztab_psh.loc[:, "modifications"] = out_mztab_psh.apply(
                lambda x: find_modification(x["opt_global_cv_MS:1000889_peptidoform_sequence"]), axis=1,
                result_type="expand"
            )

            out_mztab_psh.loc[:, "spectra_ref"] = out_mztab_psh.apply(
                lambda x: "ms_run[{}]:".format(x["ms_run"]) + x["opt_global_spectrum_reference"], axis=1,
                result_type="expand"
            )

            out_mztab_psh.loc[:, "opt_global_cv_MS:1000889_peptidoform_sequence"] = out_mztab_psh.apply(
                lambda x: AASequence.fromString(x["opt_global_cv_MS:1000889_peptidoform_sequence"]).toString(),
                axis=1,
                result_type="expand",
            )

            new_cols = [col for col in out_mztab_psh.columns if not col.startswith("opt_")] + [
                col for col in out_mztab_psh.columns if col.startswith("opt_")
            ]
            out_mztab_psh = out_mztab_psh[new_cols]
            out_mztab_psh.loc[:, 'scan_number'] = out_mztab_psh['spectra_ref'].apply(lambda x: generate_scan_number(x))
            out_mztab_psh['spectra_ref'] = out_mztab_psh['spectra_ref'].apply(lambda x: self._ms_runs[x.split(":")[0]])
            parquet_table = self.generate_psm_file(out_mztab_psh)
            
            if not pqwriter:
                # create a parquet write object giving it an output file
                pqwriter = pq.ParquetWriter(psm_output_path, parquet_table.schema)
            pqwriter.write_table(parquet_table)

            spectra_count_dict = self.__get_spectra_count(out_mztab_psh, spectra_count_dict)
            map_dict, psm_unique_keys = self.__extract_from_psm_to_pep_msg(out_mztab_psh, map_dict, psm_unique_keys)
        if pqwriter:
            pqwriter.close()
        return map_dict, spectra_count_dict
    
    def __extract_from_psm_to_pep_msg(self, psm, map_dict, psm_unique_keys):

        for key, df in psm.groupby(['opt_global_cv_MS:1000889_peptidoform_sequence', 'charge']):
            if key not in map_dict.keys():
                map_dict[key] = [None, None, None]
                psm_unique_keys.append(key)
            df = df.reset_index(drop=True)
            if pd.isna(map_dict[key][1]):
                if 'opt_global_q-value_score' in df.columns:
                    map_dict[key][0] = df.iloc[df['opt_global_q-value_score'].idxmin()]['opt_global_q-value_score']
                    map_dict[key][1] = df.iloc[df['opt_global_q-value_score'].idxmin()]['spectra_ref']
                    map_dict[key][2] = df.iloc[df['opt_global_q-value_score'].idxmin()]['scan_number']
                elif 'search_engine_score[1]' in df.columns:
                    map_dict[key][0] = df.iloc[df['search_engine_score[1]'].idxmin()]['search_engine_score[1]']
                    map_dict[key][1] = df.iloc[df['search_engine_score[1]'].idxmin()]['spectra_ref']
                    map_dict[key][2] = df.iloc[df['search_engine_score[1]'].idxmin()]['scan_number']
                else:
                    raise Exception(
                        "The psm table don't have opt_global_q-value_score or search_engine_score[1] columns")
            elif key in psm_unique_keys:
                if 'opt_global_q-value_score' in df.columns:
                    best_qvalue = df.iloc[df['opt_global_q-value_score'].idxmin()]['opt_global_q-value_score']
                    if float(map_dict[key][0]) > float(best_qvalue):
                        map_dict[key][0] = best_qvalue
                        map_dict[key][1] = df.iloc[df['opt_global_q-value_score'].idxmin()]['spectra_ref']
                        map_dict[key][2] = df.iloc[df['opt_global_q-value_score'].idxmin()]['scan_number']
                elif 'search_engine_score[1]' in df.columns:
                    best_qvalue = df.iloc[df['search_engine_score[1]'].idxmin()]['search_engine_score[1]']
                    if float(map_dict[key][0]) > float(best_qvalue):
                        map_dict[key][0] = best_qvalue
                        map_dict[key][1] = df.iloc[df['search_engine_score[1]'].idxmin()]['spectra_ref']
                        map_dict[key][2] = df.iloc[df['search_engine_score[1]'].idxmin()]['scan_number']
            if len(map_dict[key]) == 3:
                map_dict[key].append(df['start'].values[0])
                map_dict[key].append(df['end'].values[0])
                map_dict[key].append(df['unique'].values[0])
                map_dict[key].append(df['modifications'].values[0])
            if 'opt_global_Posterior_Error_Probability_score' in df.columns or 'opt_global_Posterior_Error_Probability' in df.columns:
                if len(map_dict[key]) != 7:
                    if 'opt_global_Posterior_Error_Probability_score' in df.columns:
                        probability_score = df['opt_global_Posterior_Error_Probability_score'].min()
                    else:
                        probability_score = df['opt_global_Posterior_Error_Probability'].min()
                    if float(probability_score) < map_dict[key][7]:
                        map_dict[key][7] = probability_score
                else:
                    if 'opt_global_Posterior_Error_Probability_score' in df.columns:
                        map_dict[key].append(df['opt_global_Posterior_Error_Probability_score'].min())
                    else:
                        map_dict[key].append(df['opt_global_Posterior_Error_Probability'].min())
            else:
                if len(map_dict[key]) == 7:
                    map_dict[key].append(None)
            if len(map_dict[key]) == 8:
                if "opt_global_cv_MS:1002217_decoy_peptide" in df.columns:
                    map_dict[key].append(
                        df["opt_global_cv_MS:1002217_decoy_peptide"].values[0]
                    )
                else:
                    map_dict[key].append(None)
            if len(map_dict[key]) != 11:
                if map_dict[key][1] not in df["spectra_ref"].values:
                    map_dict[key].append(df["calc_mass_to_charge"].values[0])
                    map_dict[key].append(None)
                else:
                    cals = df[
                        (df["spectra_ref"] == map_dict[key][1])
                        & (df["scan_number"] == map_dict[key][2])
                        ]["calc_mass_to_charge"].values
                    if len(cals) == 0:
                        map_dict[key].append(None)
                        map_dict[key].append(None)
                    else:
                        map_dict[key].append(cals[0])
                        map_dict[key].append(
                            df[
                                (df["spectra_ref"] == map_dict[key][1])
                                & (df["scan_number"] == map_dict[key][2])
                                ]["exp_mass_to_charge"].values[0]
                        )
            elif map_dict[key][-1] == None:
                if map_dict[key][1] in df["spectra_ref"].values:
                    cals = df[
                        (df["spectra_ref"] == map_dict[key][1])
                        & (df["scan_number"] == map_dict[key][2])
                        ]["calc_mass_to_charge"].values
                    if len(cals) != 0:
                        map_dict[key][-2] = cals[0]
                        map_dict[key][-1] = df[
                            (df["spectra_ref"] == map_dict[key][1])
                            & (df["scan_number"] == map_dict[key][2])
                            ]["exp_mass_to_charge"].values[0]

        return map_dict, psm_unique_keys

    def __get_spectra_count(self, psm, counter):
        """
        mzTab_path: mzTab file path
        psm_chunksize: the large of chunk
        return: a dict about peptide numbers
        """
        spectra_dict = (
            psm[["opt_global_cv_MS:1000889_peptidoform_sequence", "charge", "spectra_ref",]]
            .groupby(["opt_global_cv_MS:1000889_peptidoform_sequence", "charge", "spectra_ref",])
            .size()
        )
        counter.update(spectra_dict.to_dict())
        return counter

    def __check_mbr_peptide(self, reference, scan, exp_mass):
        if exp_mass is None:
            return None, None
        else:
            return reference, scan

    def generate_feature_and_psm_file(self, report_path: str, design_file: str, fasta_path: str, modifications: List,
                              pg_path: str, pr_path: str, qvalue_threshold: float, mzml_info_folder: str,
                              sdrf_path: str, output_path: str, psm_output_path,chunksize: int):

        st = time.time()
        s_data_frame, f_table = get_exp_design_dfs(design_file)
        print_estimated_time(st, "get_exp_design_dfs")

        st = time.time()
        _, index_ref = self.generete_fasta_msg(fasta_path, f_table)
        print_estimated_time(st, "generete_fasta_msg")

        st = time.time()
        if not self.uniq_masses_map:
            self.uniq_masses_map, self.unique_reference_map, self.protein_best_score_dict, self.peptide_best_score_dict = self.get_uniq_and_index_and_reference_map(
                report_path, chunksize=chunksize)
        print_estimated_time(st, "get_uniq_and_index_and_reference_map")

        st = time.time()
        if not self._ms_runs:
            self._ms_runs = self.get_ms_runs(index_ref)
            self._modifications = get_modifications(modifications[0], modifications[1])
            self._protein_map = self.get_protein_map(pg_path, self.protein_best_score_dict)
        print_estimated_time(st, "get_ms_runs, modifications and protein_map")

        st = time.time()
        map_dict = self.extract_dict_from_pep(pr_path, self.peptide_best_score_dict)
        map_dict, spectra_count_dict = self.extract_from_psm_to_pep_msg(report_path,
                                                                        qvalue_threshold,
                                                                        mzml_info_folder, self.uniq_masses_map,
                                                                        index_ref, map_dict, psm_output_path, chunksize)
        print_estimated_time(st, "extract_from_psm_to_pep_msg")

        schema = FeatureHandler()
        feature = FeatureInMemory('LFQ', schema.schema)
        feature._modifications = self._modifications
        feature._score_names = self._score_names
        pqwriter = None
        for report in self.main_report_df(report_path, qvalue_threshold, chunksize, self.uniq_masses_map):
            msstats_in = self.get_msstats_in(report, self.unique_reference_map, s_data_frame, f_table)
            msstats_in['Reference'] = msstats_in['Reference'].apply(lambda x: x.split(".")[0])
            msstats_in.loc[:, 'protein_global_qvalue'] = msstats_in['ProteinName'].apply(
                lambda x: handle_protein_map(self._protein_map, x))
            msstats_in.loc[:, 'sequence'] = msstats_in['PeptideSequence'].apply(lambda x: clean_peptidoform_sequence(x))

            no_lfq_usecols = [
                col
                for col in self._lfq_msstats_usecols
                if col not in msstats_in.columns
            ]
            for col in no_lfq_usecols:
                if col == "Channel":
                    msstats_in.loc[:, col] = "LABEL FREE SAMPLE"
                else:
                    msstats_in.loc[:, col] = None
            msstats_in = msstats_in[self._lfq_msstats_usecols]

            msstats_in = feature._map_msstats_in(msstats_in, map_dict, spectra_count_dict)
            msstats_in.loc[:, 'peptidoform'] = msstats_in[['sequence', 'modifications']].apply(
                lambda row: get_peptidoform_proforma_version_in_mztab(row['sequence'], row['modifications'],
                                                                      self._modifications), axis=1)
            msstats_in.drop(['PeptideSequence'], inplace=True, axis=1)
            msstats_in[["best_psm_reference_file_name", "best_psm_scan_number"]] = msstats_in[
                ["best_psm_reference_file_name", "best_psm_scan_number", 'exp_mass_to_charge']
            ].apply(
                lambda row: self.__check_mbr_peptide(row["best_psm_reference_file_name"], row["best_psm_scan_number"],
                                                     row['exp_mass_to_charge']),
                axis=1,
                result_type="expand",
            )
            table = feature._merge_sdrf_to_msstats_in(sdrf_path, msstats_in)
            parquet_table = feature.convert_to_parquet(table)
            if not pqwriter:
                # create a parquet write object giving it an output file
                pqwriter = pq.ParquetWriter(output_path, parquet_table.schema)
            pqwriter.write_table(parquet_table)
        if pqwriter:
            pqwriter.close()

    def generate_psm_file(self,psm):
        use_cols = [
            "sequence",
            "accession",
            "start",
            "end",
            "unique",
            "search_engine_score[1]",
            "modifications",
            "retention_time",
            "charge",
            "calc_mass_to_charge",
            "exp_mass_to_charge",
            "spectra_ref",
            'scan_number',
            "opt_global_q-value",
            "opt_global_cv_MS:1002217_decoy_peptide"
        ]
        psm = psm[use_cols].copy()
        psm.loc[:, 'protein_global_qvalue'] = psm['accession'].apply(
            lambda x: handle_protein_map(self._protein_map, x))
        psm.loc[:, 'peptidoform'] = psm[['sequence', 'modifications']].apply(
            lambda row: get_peptidoform_proforma_version_in_mztab(row['sequence'], row['modifications'],
                                                                    self._modifications), axis=1)
        psm.loc[:, 'id_scores'] = psm['search_engine_score[1]'].apply(lambda x: [
            f"protein-level q-value: {x}",
            f"Posterior error probability: ",
        ])
        psm.drop(["search_engine_score[1]"], inplace=True, axis=1)

        psm.loc[:, "posterior_error_probability"] = None
        psm.loc[:, "consensus_support"] = None

        psm.rename(columns=self.map_psm, inplace=True)
        parquet_table = self.convert_psm_format_to_parquet(psm)

        return parquet_table

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
        res['protein_start_positions'] = res['protein_start_positions'].apply(
            self.__split_start_or_end).to_list()
        res['protein_end_positions'] = res['protein_end_positions'].apply(
            self.__split_start_or_end).to_list()
        res['protein_global_qvalue'] = res['protein_global_qvalue'].astype(float)
        res['unique'] = res['unique'].map(lambda x: pd.NA if pd.isna(x) else int(x)).astype('Int32')
        res['modifications'] = res['modifications'].apply(lambda x: feature._generate_modification_list(x))
        res['retention_time'] = res['retention_time'].astype(float)
        res['charge'] = res['charge'].map(lambda x: pd.NA if pd.isna(x) else int(x)).astype('Int32')
        res['exp_mass_to_charge'] = res['exp_mass_to_charge'].astype(float)
        res['calc_mass_to_charge'] = res['calc_mass_to_charge'].astype(float)
        res['posterior_error_probability'] = res['posterior_error_probability'].astype(float)

        res['is_decoy'] = res['is_decoy'].map(lambda x: pd.NA if pd.isna(x) else int(x)).astype('Int32')

        res.loc[:, "num_peaks"] = None
        res.loc[:, "mz_array"] = None
        res.loc[:, "intensity_array"] = None

        res.loc[:, "gene_accessions"] = None
        res.loc[:, "gene_names"] = None

        return pa.Table.from_pandas(res, schema=psm._create_schema())

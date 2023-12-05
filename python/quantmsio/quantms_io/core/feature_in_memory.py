import codecs
from quantms_io.core.mztab import fetch_modifications_from_mztab_line
from quantms_io.utils.pride_utils import clean_peptidoform_sequence, get_petidoform_msstats_notation,get_quantmsio_modifications,get_peptidoform_proforma_version_in_mztab,generate_scan_number

import numpy as np
import pandas as pd
import os
import pyarrow as pa
import pyarrow.parquet as pq
from quantms_io.utils.constants import ITRAQ_CHANNEL, TMT_CHANNELS
import swifter

def get_modifications(mztab_path):
    """
    mzTab_path: mzTab file path
    return: a dict about modifications
    """
    if os.stat(mztab_path).st_size == 0:
        raise ValueError("File is empty")
    f = codecs.open(mztab_path, "r", "utf-8")
    line = f.readline()
    mod_dict = {}
    while line.split("\t")[0] == "MTD":
        if "_mod[" in line:
            mod_dict = fetch_modifications_from_mztab_line(line, mod_dict)
        line = f.readline()
    f.close()
    return mod_dict



class FeatureInMemory:
    def __init__(self, experiment_type, schema):
        self.mzml_directory = None
        self.experiment_type = experiment_type
        self.schema = schema
        # mzTab file
        self._mzTab_file = None
        # psm pos
        self._psm_pos = None
        # psm len
        self._psm_len = None
        # pep pos
        self._pep_pos = None
        # pep len
        self._pep_len = None
        # prt pos
        self._prt_pos = None
        # prt len
        self._prt_len = None
        # load psms columns
        self._psms_columns = None
        # load pep columns
        self._pep_columns = None
        # ms_runs
        self._ms_runs = None
        # msstats_in usecols
        self._tmt_msstats_usecols = [
            "ProteinName",
            "PeptideSequence",
            "Channel",
            "Run",
            "BioReplicate",
            "Intensity",
            "FragmentIon",
            "IsotopeLabelType",
            "Charge",
            "Reference",
            "protein_global_qvalue",
            "sequence",
        ]
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
        # score names
        self._score_names = None  # dict
        # map dict
        self._map_lfq = {
            'ProteinName': 'protein_accessions',
            'Reference': 'reference_file_name',
            'Run': 'run',
            'BioReplicate': 'biological_replicate',
            'Intensity': 'intensity',
            'FragmentIon': 'fragment_ion',
            'IsotopeLabelType': 'isotope_label_type',
            'PrecursorCharge': 'charge',
            'Channel': 'channel',
            'source name': 'sample_accession',
            'comment[fraction identifier]': 'fraction',
        }
        self._map_tmt = {
            'ProteinName': 'protein_accessions',
            'Reference': 'reference_file_name',
            'Run': 'run',
            'BioReplicate': 'biological_replicate',
            'Intensity': 'intensity',
            'FragmentIon': 'fragment_ion',
            'IsotopeLabelType': 'isotope_label_type',
            'RetentionTime': 'retention_time',
            'Charge': 'charge',
            'Channel': 'channel',
            'source name': 'sample_accession',
            'comment[fraction identifier]': 'fraction',
        }

    def __set_table_config(self, header, length, pos, file_name):
        """
        set table config.
        header: table tag #PSH PRH PEH
        length: table length
        pos: table pos
        file_name: mzTab file name
        """
        self._mzTab_file = file_name
        if header == "PSH":
            self._psm_pos = pos
            self._psm_len = length
        elif header == "PEH":
            self._pep_pos = pos
            self._pep_len = length
        else:
            self._prt_pos = pos
            self._prt_len = length

    # second load
    def __load_second(self, fle, header, **kwargs):
        f = open(fle)
        if header == "PSH":
            f.seek(self._psm_pos)
            return pd.read_csv(f, nrows=self._psm_len, **kwargs)
        elif header == "PEH":
            f.seek(self._pep_pos)
            return pd.read_csv(f, nrows=self._pep_len, **kwargs)
        else:
            f.seek(self._prt_pos)
            return pd.read_csv(f, nrows=self._prt_len, **kwargs)

    # extract pep columns
    def __extract_pep_columns(self, file):
        if os.stat(file).st_size == 0:
            raise ValueError("File is empty")
        f = open(file)
        line = f.readline()
        while not line.startswith("PEH"):
            line = f.readline()
        self._pep_columns = line.split("\n")[0].split("\t")

    # extract csv len
    def _extract_len(self, fle, header):
        map_tag = {"PSH": "PSM", "PEH": "PEP", "PRH": "PRT"}
        if os.stat(fle).st_size == 0:
            raise ValueError("File is empty")
        f = open(fle)
        pos = 0
        line = f.readline()
        while not line.startswith(header):
            pos = f.tell()
            line = f.readline()

        if header == "PSH":
            self._psms_columns = line.split("\n")[0].split("\t")

        line = f.readline()
        fle_len = 0
        while line.startswith(map_tag[header]):
            fle_len += 1
            line = f.readline()
        f.close()
        return fle_len, pos

    def skip_and_load_csv(self, fle, header, **kwargs):
        """
        file: mzTab file
        :param fle: mzTab file
        :param header table tag #PSH PRH PEH
        """
        if self._mzTab_file == fle and self._psm_pos is not None and header == "PSH":
            return self.__load_second(fle, header, **kwargs)
        if self._mzTab_file == fle and self._pep_pos is not None and header == "PEH":
            return self.__load_second(fle, header, **kwargs)
        if self._mzTab_file == fle and self._prt_pos is not None and header == "PRH":
            return self.__load_second(fle, header, **kwargs)
        fle_len, pos = self._extract_len(fle, header)
        if os.stat(fle).st_size == 0:
            raise ValueError("File is empty")
        f = open(fle)
        f.seek(pos)
        self.__set_table_config(header, fle_len, pos, fle)
        return pd.read_csv(f, nrows=fle_len, **kwargs)

    def __get_spectra_count(self, mztab_path, psm_chunksize):
        """
        mzTab_path: mzTab file path
        psm_chunksize: the large of chunk
        return: a dict about piptie numbers
        """
        from collections import Counter

        counter = Counter()
        psms = self.skip_and_load_csv(
            mztab_path, "PSH", sep="\t", chunksize=psm_chunksize
        )
        for psm in psms:
            psm["spectra_ref"] = psm["spectra_ref"].swifter.apply(
                lambda x: self._ms_runs[x.split(":")[0]]
            )
            if "opt_global_cv_MS:1000889_peptidoform_sequence" not in psm.columns:
                psm.loc[:, 'opt_global_cv_MS:1000889_peptidoform_sequence'] = psm[['modifications', 'sequence']].swifter.apply(
                    lambda row: get_petidoform_msstats_notation(row['sequence'], row['modifications'], self._modifications),
                    axis=1)
            spectra_dict = (
                psm[
                    [
                        "opt_global_cv_MS:1000889_peptidoform_sequence",
                        "charge",
                        "spectra_ref",
                    ]
                ]
                .groupby(
                    [
                        "opt_global_cv_MS:1000889_peptidoform_sequence",
                        "charge",
                        "spectra_ref",
                    ]
                )
                .size()
            )
            counter.update(spectra_dict.to_dict())
        return counter

    def _get_protein_map(self, mztab_path):
        """
        return: a dict about protein score
        """
        prt = self.skip_and_load_csv(
            mztab_path,
            "PRH",
            sep="\t",
            usecols=["ambiguity_members", "best_search_engine_score[1]"],
        )
        prt_score = prt.groupby("ambiguity_members").min()
        protein_map = prt_score.to_dict()["best_search_engine_score[1]"]
        return protein_map

    @staticmethod
    def _get_score_names(fle):
        """
        return: a dict about search engine
        """
        if os.stat(fle).st_size == 0:
            raise ValueError("File is empty")
        f = codecs.open(fle, "r", "utf-8")
        line = f.readline()
        score_names = {}
        while line.split("\t")[0] == "MTD":
            if "search_engine_score" in line:
                score_values = line.replace("[", "").replace("]", "").split(",")
                score_name = score_values[2].strip()
                if ":" in score_name:
                    score_name = "'{}'".format(score_name)
                if "peptide" in line.split("\t")[1]:
                    score_names["peptide_score"] = score_name
                elif "protein" in line.split("\t")[1]:
                    score_names["protein_score"] = score_name
                else:
                    score_names["psm_score"] = score_name

            line = f.readline()
        f.close()
        return score_names

    @staticmethod
    def __handle_protein_map(protein_map, key):
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

    def _extract_from_pep(self, mztab_path):
        """
        return: dict about pep_msg
        """
        self.__extract_pep_columns(mztab_path)
        pep_usecols = [
            "opt_global_cv_MS:1000889_peptidoform_sequence",
            "charge",
            "best_search_engine_score[1]",
            "spectra_ref",
        ]
        live_cols = [col for col in pep_usecols if col in self._pep_columns]
        not_cols = [col for col in pep_usecols if col not in live_cols]
        if "opt_global_cv_MS:1000889_peptidoform_sequence" in not_cols:
            if "sequence" in self._pep_columns and "modifications" in self._pep_columns:
                live_cols.append("sequence")
                live_cols.append("modifications")
            else:
                raise Exception(
                    "The peptide table don't have opt_global_cv_MS:1000889_peptidoform_sequence columns"
                )
        if "charge" in not_cols or "best_search_engine_score[1]" in not_cols:
            raise Exception(
                "The peptide table don't have best_search_engine_score[1] or charge columns"
            )

        pep = self.skip_and_load_csv(mztab_path, "PEH", sep="\t", usecols=live_cols)

        # check opt_global_cv_MS:1000889_peptidoform_sequence
        if "opt_global_cv_MS:1000889_peptidoform_sequence" not in pep.columns:
            modifications = get_modifications(mztab_path)
            pep.loc[:, "opt_global_cv_MS:1000889_peptidoform_sequence"] = pep[
                ["modifications", "sequence"]
            ].swifter.apply(
                lambda row: get_petidoform_msstats_notation(
                    row["sequence"], row["modifications"], modifications
                ),
                axis=1,
            )

        # check spectra_ref
        if "spectra_ref" not in pep.columns:
            pep.loc[:, "scan_number"] = None
            pep.loc[:, "spectra_ref"] = None
        else:
            pep.loc[:, "scan_number"] = pep["spectra_ref"].swifter.apply(
                lambda x: generate_scan_number(x)
            )
            pep["spectra_ref"] = pep["spectra_ref"].swifter.apply(
                lambda x: self._ms_runs[x.split(":")[0]]
            )
        pep_msg = pep.iloc[
            pep.groupby(
                ["opt_global_cv_MS:1000889_peptidoform_sequence", "charge"]
            ).swifter.apply(lambda row: row["best_search_engine_score[1]"].idxmin())
        ]
        pep_msg = pep_msg.set_index(
            ["opt_global_cv_MS:1000889_peptidoform_sequence", "charge"]
        )

        pep_msg.loc[:, "pep_msg"] = pep_msg[
            ["best_search_engine_score[1]", "spectra_ref", "scan_number"]
        ].swifter.apply(
            lambda row: [
                row["best_search_engine_score[1]"],
                row["spectra_ref"],
                row["scan_number"],
            ],
            axis=1,
        )

        map_dict = pep_msg.to_dict()["pep_msg"]
        return map_dict

    def _extract_from_psm_to_pep_msg(self, mztab_path, map_dict):
        """
        return dict about pep and psm msg
        """
        psms = self.skip_and_load_csv(mztab_path, 'PSH', sep='\t',dtype={'start':str,'end':str},chunksize=1000000)
        self._modifications = get_modifications(mztab_path)
        psm_unique_keys =[]
        for psm in psms:
            if 'opt_global_cv_MS:1000889_peptidoform_sequence' not in psm.columns:
                psm.loc[:, 'opt_global_cv_MS:1000889_peptidoform_sequence'] = psm[['modifications', 'sequence']].swifter.apply(
                    lambda row: get_petidoform_msstats_notation(row['sequence'], row['modifications'], self._modifications),
                    axis=1)
            for key, df in psm.groupby(['opt_global_cv_MS:1000889_peptidoform_sequence', 'charge']):
                if key not in map_dict.keys():
                    map_dict[key] = [None, None, None]
                    psm_unique_keys.append(key)
                df = df.reset_index(drop=True)
                df.loc[:, 'scan_number'] = df['spectra_ref'].swifter.apply(lambda x: generate_scan_number(x))
                df['spectra_ref'] = df['spectra_ref'].swifter.apply(lambda x: self._ms_runs[x.split(":")[0]])
                if pd.isna(map_dict[key][1]):
                    if 'opt_global_q-value_score' in df.columns:
                        temp_df = df.iloc[df['opt_global_q-value_score'].idxmin()]
                        map_dict[key][0] = temp_df['opt_global_q-value_score']
                        map_dict[key][1] = temp_df['spectra_ref']
                        map_dict[key][2] = temp_df['scan_number']
                    elif 'search_engine_score[1]' in df.columns:
                        temp_df = df.iloc[df['search_engine_score[1]'].idxmin()]
                        map_dict[key][0] = temp_df['search_engine_score[1]']
                        map_dict[key][1] = temp_df['spectra_ref']
                        map_dict[key][2] = temp_df['scan_number']
                    else:
                        raise Exception(
                            "The psm table don't have opt_global_q-value_score or search_engine_score[1] columns")
                elif key in psm_unique_keys:
                    if 'opt_global_q-value_score' in df.columns:
                        temp_df = df.iloc[df['opt_global_q-value_score'].idxmin()]
                        best_qvalue = temp_df['opt_global_q-value_score']
                        if float(map_dict[key][0]) > float(best_qvalue):
                            map_dict[key][0] = best_qvalue
                            map_dict[key][1] = temp_df['spectra_ref']
                            map_dict[key][2] = temp_df['scan_number']
                    elif 'search_engine_score[1]' in df.columns:
                        temp_df = df.iloc[df['search_engine_score[1]'].idxmin()]
                        best_qvalue = temp_df['search_engine_score[1]']
                        if float(map_dict[key][0]) > float(best_qvalue):
                            map_dict[key][0] = best_qvalue
                            map_dict[key][1] = temp_df['spectra_ref']
                            map_dict[key][2] = temp_df['scan_number']
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
                            ]
                        
                        if len(cals) == 0:
                            map_dict[key].append(None)
                            map_dict[key].append(None)
                        else:
                            map_dict[key].append(cals["calc_mass_to_charge"].values[0])
                            map_dict[key].append(
                                cals["exp_mass_to_charge"].values[0]
                            )
                elif map_dict[key][-1] == None:
                    if map_dict[key][1] in df["spectra_ref"].values:
                        cals = df[
                            (df["spectra_ref"] == map_dict[key][1])
                            & (df["scan_number"] == map_dict[key][2])
                            ]
                        if len(cals) != 0:
                            map_dict[key][-2] = cals["calc_mass_to_charge"].values[0]
                            map_dict[key][-1] = cals["exp_mass_to_charge"].values[0]
                        
        return map_dict

    def _extract_psm_pep_msg(self, mztab_path):
        """
        mzTab_path: mzTab file path
        return: dict about pep and psm msg
        """
        # load ms_runs
        self._ms_runs = self.extract_ms_runs(mztab_path)
        self._score_names = self._get_score_names(mztab_path)
        map_dict = self._extract_from_pep(mztab_path)
        map_dict = self._extract_from_psm_to_pep_msg(mztab_path, map_dict)

        return map_dict

    def _extract_rt_from_consensus_xml(self,intensity_map,msstats_in):
        if "RetentionTime" not in msstats_in.columns:
            if self.experiment_type != 'LFQ':
                msstats_in['retention_time'] = msstats_in[['peptidoform','Charge','Reference','Channel','Intensity']].swifter.apply(
                    lambda row: self.__map_rt_or_exp_mass(row[:-1],intensity_map,row[-1:].values[0],label='rt'), axis=1
                )
            else:
                msstats_in['retention_time'] = msstats_in[['peptidoform','PrecursorCharge','Reference','Intensity']].swifter.apply(
                    lambda row: self.__map_rt_or_exp_mass(row[:-1],intensity_map,row[-1:].values[0],label='rt'), axis=1
                )
        if self.experiment_type != 'LFQ':
            msstats_in['exp_mass_to_charge'] = msstats_in[['peptidoform','Charge','Reference','Channel','Intensity','exp_mass_to_charge']].swifter.apply(
                    lambda row: self.__map_rt_or_exp_mass(row[:-2],intensity_map,row[-2:-1].values[0],label='exp',exp_mass=row[-1:].values[0]), axis=1
            )
        else:
            msstats_in['exp_mass_to_charge'] = msstats_in[['peptidoform','PrecursorCharge','Reference','Intensity','exp_mass_to_charge']].swifter.apply(
                    lambda row: self.__map_rt_or_exp_mass(row[:-2],intensity_map,row[-2:-1].values[0],label='exp',exp_mass=row[-1:].values[0]), axis=1
            )
        return msstats_in

    def __map_rt_or_exp_mass(self, row, intensity_map, intensity, label=None, exp_mass=None):
        row = list(map(str, row.tolist()))
        key = ":_:".join(row)
        if key in intensity_map:
            if abs(intensity_map[key]["intensity"] - np.float64(intensity)) < 0.1:
                if label == 'rt':
                    return intensity_map[key]["rt"]
                else:
                    return intensity_map[key]["mz"]
            else:
                if label =='exp':
                    return exp_mass
                else:
                    return None
        else:
            if label == 'exp':
                return exp_mass
            else:
                return None

    def __check_mbr_peptide(self,reference,scan,exp_mass):
        if exp_mass is None:
            return None,None
        else:
            return reference,scan

    def _map_msstats_in(self, msstats_in, map_dict, spectra_count_dict):
        """
        map key: PeptideSequence, charge
        :param msstats_in: msstats_in dataframe
        :param map_dict: dict about pep and psm msg
        :param spectra_count_dict: dict about piptie numbers
        :return: msstats_in dataframe
        """
        if self.experiment_type == "LFQ":
            msstats_in.loc[:, "spectral_count"] = msstats_in[
                ["PeptideSequence", "PrecursorCharge", "Reference"]
            ].swifter.apply(
                lambda row: spectra_count_dict[
                    (row["PeptideSequence"], row["PrecursorCharge"], row["Reference"])
                ],
                axis=1,
            )
        else:
            msstats_in.loc[:, "spectral_count"] = msstats_in[
                ["PeptideSequence", "Charge", "Reference"]
            ].swifter.apply(
                lambda row: spectra_count_dict[
                    (row["PeptideSequence"], row["Charge"], row["Reference"])
                ],
                axis=1,
            )
        map_features = [
            "global_qvalue",
            "best_psm_reference_file_name",
            "best_psm_scan_number",
            "protein_start_positions",
            "protein_end_positions",
            "unique",
            "modifications",
            "posterior_error_probability",
            "is_decoy",
            "calc_mass_to_charge",
            "exp_mass_to_charge",
        ]
        for i, feature in enumerate(map_features):
            if self.experiment_type == "LFQ":
                msstats_in.loc[:, feature] = msstats_in[
                    ["PeptideSequence", "PrecursorCharge"]
                ].swifter.apply(
                    lambda row: map_dict[
                        (row["PeptideSequence"], row["PrecursorCharge"])
                    ][i],
                    axis=1,
                )
            else:
                msstats_in.loc[:, feature] = msstats_in[
                    ["PeptideSequence", "Charge"]
                ].swifter.apply(
                    lambda row: map_dict[(row["PeptideSequence"], row["Charge"])][i],
                    axis=1,
                )
        peptide_score_name = self._score_names["peptide_score"]
        msstats_in["id_scores"] = (
            peptide_score_name + ":" + msstats_in["global_qvalue"].astype(str) + ","
            + "Best PSM PEP:" + msstats_in["posterior_error_probability"].astype(str)
        )
        return msstats_in

    def merge_mztab_and_sdrf_to_msstats_in(
        self,
        mztab_path,
        msstats_path,
        sdrf_path,
        output_path,
        msstats_chunksize=1000000,
        intensity_map = None
    ):
        """
        mzTab_path: mzTab file path
        msstats_path: msstats_in file path
        sdrf_path: sdrf file path
        output_path: output path of parquet file
        msstats_chunksize: the large of msstats chunk
        """
        protein_map = self._get_protein_map(mztab_path)
        map_dict = self._extract_psm_pep_msg(mztab_path)
        msstats_ins = pd.read_csv(msstats_path, chunksize=msstats_chunksize)
        spectra_count_dict = self.__get_spectra_count(mztab_path, 500000)
        pqwriter = None
        for msstats_in in msstats_ins:
            msstats_in['Reference'] = msstats_in['Reference'].swifter.apply(
                lambda x: x.split(".")[0])
            msstats_in.loc[:, 'protein_global_qvalue'] = msstats_in['ProteinName'].swifter.apply(
                lambda x: self.__handle_protein_map(protein_map, x))
            msstats_in.loc[:, 'sequence'] = (msstats_in['PeptideSequence']
                                             .swifter.apply(lambda x: clean_peptidoform_sequence(x)))
            
            if self.experiment_type != 'LFQ':
                no_tmt_usecols = [
                    col
                    for col in self._tmt_msstats_usecols
                    if col not in msstats_in.columns
                ]
                for col in no_tmt_usecols:
                    if "IsotopeLabelType" == col:
                        msstats_in.loc[:, col] = "L"
                    else:
                        msstats_in.loc[:, col] = None
                if "TMT" in self.experiment_type:
                    msstats_in["Channel"] = msstats_in["Channel"].swifter.apply(
                        lambda row: TMT_CHANNELS[self.experiment_type][row - 1]
                    )
                else:
                    msstats_in["Channel"] = msstats_in["Channel"].swifter.apply(
                        lambda row: ITRAQ_CHANNEL[self.experiment_type][row - 1]
                    )
                if "RetentionTime" in msstats_in.columns:
                    self._tmt_msstats_usecols.append("RetentionTime")
                msstats_in = msstats_in[self._tmt_msstats_usecols]
                self._tmt_msstats_usecols.remove('RetentionTime')
                msstats_in = self._map_msstats_in(msstats_in, map_dict, spectra_count_dict)
                msstats_in.loc[:,'peptidoform'] = msstats_in[['sequence','modifications']].swifter.apply(lambda row: get_peptidoform_proforma_version_in_mztab(row['sequence'],row['modifications'],self._modifications),axis=1)
                msstats_in.drop(['PeptideSequence'],inplace=True, axis=1)
                msstats_in[["best_psm_reference_file_name", "best_psm_scan_number"]] = msstats_in[
                ["best_psm_reference_file_name", "best_psm_scan_number",'exp_mass_to_charge']
                ].swifter.apply(
                lambda row: self.__check_mbr_peptide(row["best_psm_reference_file_name"],row["best_psm_scan_number"],row['exp_mass_to_charge']),
                axis=1,
                result_type="expand",
                )
                if intensity_map is not None and len(intensity_map)!=0:
                    msstats_in = self._extract_rt_from_consensus_xml(intensity_map,msstats_in)
                table = self._merge_sdrf_to_msstats_in(sdrf_path, msstats_in)
                parquet_table = self.convert_to_parquet(table)
                if not pqwriter:
                    # create a parquet write object giving it an output file
                    pqwriter = pq.ParquetWriter(output_path, parquet_table.schema)
                pqwriter.write_table(parquet_table)
            else:
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
                msstats_in = self._map_msstats_in(msstats_in, map_dict, spectra_count_dict)
                msstats_in.loc[:,'peptidoform'] = msstats_in[['sequence','modifications']].swifter.apply(lambda row: get_peptidoform_proforma_version_in_mztab(row['sequence'],row['modifications'],self._modifications),axis=1)
                msstats_in.drop(['PeptideSequence'],inplace=True, axis=1)
                msstats_in[["best_psm_reference_file_name", "best_psm_scan_number"]] = msstats_in[
                ["best_psm_reference_file_name", "best_psm_scan_number",'exp_mass_to_charge']
                ].swifter.apply(
                lambda row: self.__check_mbr_peptide(row["best_psm_reference_file_name"],row["best_psm_scan_number"],row['exp_mass_to_charge']),
                axis=1,
                result_type="expand",
                )
                if intensity_map is not None and len(intensity_map)!=0:
                    msstats_in = self._extract_rt_from_consensus_xml(intensity_map,msstats_in)
                table = self._merge_sdrf_to_msstats_in(sdrf_path, msstats_in)
                parquet_table = self.convert_to_parquet(table)
                if not pqwriter:
                    # create a parquet write object giving it an output file
                    pqwriter = pq.ParquetWriter(output_path, parquet_table.schema)
                pqwriter.write_table(parquet_table)
        if pqwriter:
            pqwriter.close()

    def _merge_sdrf_to_msstats_in(self, sdrf_path, msstats_in):
        """
        sdrf_path: sdrf file path
        msstats_in: msstats_in dataframe
        output_path: output file path(csv)
        header: whether to write header
        """
        sdrf = pd.read_csv(sdrf_path, sep="\t")
        factor = "".join(filter(lambda x: x.startswith("factor"), sdrf.columns))
        sdrf = sdrf[
            [
                "comment[data file]",
                "source name",
                factor,
                "comment[fraction identifier]",
                "comment[label]",
            ]
        ]
        sdrf["comment[data file]"] = sdrf["comment[data file]"].swifter.apply(
            lambda x: x.split(".")[0]
        )
        if self.experiment_type != "LFQ":
            res = pd.merge(
                msstats_in,
                sdrf,
                left_on=["Reference", "Channel"],
                right_on=["comment[data file]", "comment[label]"],
                how="left",
            )
            res.drop(["comment[data file]", "comment[label]"], axis=1, inplace=True)
            res.rename(columns=self._map_tmt, inplace=True)
            res.rename(columns={factor: 'condition'},inplace=True)
            return res
        else:
            res = pd.merge(
                msstats_in,
                sdrf,
                left_on=["Reference"],
                right_on=["comment[data file]"],
                how="left",
            )
            res.drop(["comment[data file]", "comment[label]"], axis=1, inplace=True)
            res.rename(columns=self._map_lfq, inplace=True)
            res.rename(columns={factor: 'condition'},inplace=True)
            return res

    # extract ms runs
    @staticmethod
    def extract_ms_runs(fle):
        """
        fle: mzTab file path
        """
        if os.stat(fle).st_size == 0:
            raise ValueError("File is empty")
        f = codecs.open(fle, "r", "utf-8")
        line = f.readline()
        ms_runs = {}
        while line.split("\t")[0] == "MTD":
            if line.split("\t")[1].split("-")[-1] == "location":
                ms_runs[line.split("\t")[1].split("-")[0]] = (
                    line.split("\t")[2].split("//")[-1].split(".")[0]
                )
            line = f.readline()
        f.close()
        return ms_runs

    @staticmethod
    def __split_start_or_end(value):
        """
        split start or end
        :param value: start or end
        """
        if pd.isna(value) or value=='null':
            return pd.NA
        elif "," in str(value):
            return list(map(int, value.split(",")))
        elif value is np.nan:
            return None
        else:
            return [int(value)]
    
    def _generate_modification_list(self, modification_str:str):

        if pd.isna(modification_str):
            return pd.NA
        modifications = get_quantmsio_modifications(modification_str,self._modifications)
        modifications_string = ""
        for key, value in modifications.items():
            modifications_string += "|".join(map(str, value["position"]))
            modifications_string = modifications_string + "-" + value["unimod_accession"] + ","
        modifications_string = modifications_string[:-1]  # Remove last comma
        modification_list =  modifications_string.split(",")

        return modification_list

    def convert_to_parquet(self, res):
        """
        res: msstats_in dataframe
        return: parquet table
        """
        if res['id_scores'].dtype == 'str':
            res['id_scores'] = res['id_scores'].str.split(',')
        res['sequence'] = res['sequence'].astype(str)
        res['protein_accessions'] = res['protein_accessions'].str.split(";")
        res['protein_start_positions'] = res['protein_start_positions'].swifter.apply(
            self.__split_start_or_end).to_list()
        res['protein_end_positions'] = res['protein_end_positions'].swifter.apply(
            self.__split_start_or_end).to_list()
        res['protein_global_qvalue'] = res['protein_global_qvalue'].astype(float)
        res['unique'] = res['unique'].map(lambda x: pd.NA if pd.isna(x) else int(x)).astype('Int32')
        res['modifications'] = res['modifications'].swifter.apply(lambda x: self._generate_modification_list(x))
        res['charge'] = res['charge'].map(lambda x: pd.NA if pd.isna(x) else int(x)).astype('Int32')
        res['exp_mass_to_charge'] = res['exp_mass_to_charge'].astype(float)
        res['calc_mass_to_charge'] = res['calc_mass_to_charge'].astype(float)
        res['posterior_error_probability'] = res['posterior_error_probability'].astype(float)
        res['global_qvalue'] = res['global_qvalue'].astype(float)
        res['is_decoy'] = res['is_decoy'].map(lambda x: pd.NA if pd.isna(x) else int(x)).astype('Int32')
        res['intensity'] = res['intensity'].astype(float)
        res['spectral_count'] = res['spectral_count'].astype(int)
        res['fraction'] = res['fraction'].astype(int).astype(str)
        res['biological_replicate'] = res['biological_replicate'].astype(str)
        res['fragment_ion'] = res['fragment_ion'].astype(str)
        res['run'] = res['run'].astype(str)
        res['best_psm_scan_number'] = res['best_psm_scan_number'].astype(str)

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

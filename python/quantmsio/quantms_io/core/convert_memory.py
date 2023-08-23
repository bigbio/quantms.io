import pyarrow as pa
import pandas as pd
import numpy as np
import os
import codecs
import re

'''
example
Convert = FeatureConvertor('lfq',schema)
Convert.merge_psm_to_pep("lfq2\PXD002854-serum.sdrf_openms_design_openms.mzTab",'res1.txt')
Convert.merge_pep_to_msstats_in("res1.txt","lfq2\PXD002854-serum.sdrf_openms_design_msstats_in.csv","res2.csv")
Convert.merge_sdrf_to_msstats_in("lfq2\PXD002854-serum.sdrf.tsv","res2.csv","result_lfq.csv")
parquet = Convert.convert_to_parquet("result_lfq.csv")
'''


tmt = {
    'TMT10':  ['TMT126', 'TMT127C', 'TMT127N', 'TMT128C', 'TMT128N', 'TMT129C', 'TMT129N', 'TMT130C', 'TMT130N', 'TMT131'],
    'TMT11': ["TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N", "TMT130C", "TMT131N", "TMT131C"],
    'TMT16': ["TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N", "TMT130C", "TMT131N", "TMT131C", "TMT132N", "TMT132C", "TMT133N", "TMT133C", "TMT134N"],
    'TMT6': ["TMT126", "TMT127", "TMT128", "TMT129", "TMT130", "TMT131"]
}


class FeatureConvertor():
    def __init__(self, experiment_type, schema):
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
        # psm usecols
        self._tmt_usecols = ['sequence', 'accession', 'unique', 'search_engine_score[1]', 'modifications', 'retention_time', 'charge', 'exp_mass_to_charge', 'calc_mass_to_charge', 'start', 'end',
                             'opt_global_Posterior_Error_Probability_score', 'opt_global_cv_MS:1000889_peptidoform_sequence', 'opt_global_cv_MS:1002217_decoy_peptide', 'opt_global_consensus_support', 'spectra_ref']
        self._lfq_usecols = ['sequence', 'accession', 'unique', 'modifications', 'retention_time', 'charge', 'exp_mass_to_charge', 'calc_mass_to_charge', 'start', 'end',
                             'opt_global_Posterior_Error_Probability_score', 'opt_global_cv_MS:1000889_peptidoform_sequence', 'opt_global_cv_MS:1002217_decoy_peptide', 'opt_global_q-value', 'spectra_ref']
        # msstats_in usecols
        self._tmt_msstats_usecols = ['ProteinName', 'PeptideSequence', 'Channel', 'Run',
                                     'BioReplicate', 'Intensity', 'FragmentIon', 'IsotopeLabelType', 'RetentionTime', 'Charge']
        self._lfq_msstats_usecols = ['ProteinName', 'PeptideSequence', 'Reference', 'Run',
                                     'BioReplicate', 'Intensity', 'FragmentIon', 'IsotopeLabelType', 'PrecursorCharge']
        # map dict
        self._map_lfq = {
            'ProteinName': 'protein_accessions',
            'PeptideSequence': 'peptidoform',
            'Reference': 'reference_file_name',
            'Run': 'run',
            'BioReplicate': 'biological_replicate',
            'Intensity': 'intensity',
            'FragmentIon': 'fragment_ion',
            'IsotopeLabelType': 'isotope_label_type',
            'PrecursorCharge': 'charge',
            'Channel': 'channel',
            'best_search_engine_score[1]': 'best_id_score',
            'search_engine_score[1]': '',
            'start': 'protein_start_positions',
            'end': 'protein_end_positions',
            'opt_global_Posterior_Error_Probability_score': 'posterior_error_probability',
            'opt_global_q-value': 'global_qvalue',
            'opt_global_cv_MS:1002217_decoy_peptide': 'is_decoy',
            'source name': 'sample_accession',
            'comment[fraction identifier]': 'fraction',
            'factor value[organism part]': 'condition'
        }
        self._map_tmt = {
            'ProteinName': 'protein_accessions',
            'PeptideSequence': 'peptidoform',
            'spectra_ref':  'reference_file_name',
            'Run': 'run',
            'BioReplicate': 'biological_replicate',
            'Intensity': 'intensity',
            'FragmentIon': 'fragment_ion',
            'IsotopeLabelType': 'isotope_label_type',
            'RetentionTime': 'retention_time',
            'Charge': 'charge',
            'Channel': 'channel',
            'best_search_engine_score[1]': 'best_id_score',
            'search_engine_score[1]': 'global_qvalue',
            'start': 'protein_start_positions',
            'end': 'protein_end_positions',
            'opt_global_Posterior_Error_Probability_score': 'posterior_error_probability',
            'opt_global_consensus_support': 'consensus_support',
            'opt_global_cv_MS:1002217_decoy_peptide': 'is_decoy',
            'source name': 'sample_accession',
            'comment[fraction identifier]': 'fraction',
            'factor value[organism part]': 'condition'
        }

    def __set_table_config(self, header, length, pos, file_name):
        self._mzTab_file = file_name
        if header == 'PSH':
            self._psm_pos = pos
            self._psm_len = length
        elif header == 'PEH':
            self._pep_pos = pos
            self._pep_len = length
        else:
            self._prt_pos = pos
            self._prt_len = length

    # second load
    def __load_second(self, fle, header, **kwargs):
        f = open(fle)
        if header == 'PSH':
            f.seek(self._psm_pos)
            return pd.read_csv(f, nrows=self._psm_len, **kwargs)
        elif header == 'PEH':
            f.seek(self._pep_pos)
            return pd.read_csv(f, nrows=self._pep_len, **kwargs)
        else:
            f.seek(self._prt_pos)
            return pd.read_csv(f, nrows=self._prt_len, **kwargs)

        # extract csv len
    def __extract_len(self, fle, header):
        map_tag = {
            "PSH": 'PSM',
            "PEH": 'PEP',
            'PRH': 'PRT'
        }
        if os.stat(fle).st_size == 0:
            raise ValueError("File is empty")
        f = codecs.open(fle, 'r', 'utf-8')
        line = f.readline()
        while line.split("\t")[0] != header:
            line = f.readline()
        line = f.readline()
        fle_len = 0
        while line.split("\t")[0] == map_tag[header]:
            fle_len += 1
            line = f.readline()
        f.close()
        return fle_len

    def __get_spectra_count(self, mzTab_path, psm_chunksize):
        from collections import Counter
        counter = Counter()
        psms = self.skip_and_load_csv(
            mzTab_path, 'PSH', sep='\t', chunksize=psm_chunksize)
        for psm in psms:
            spectra_dict = psm[['accession', 'opt_global_cv_MS:1000889_peptidoform_sequence', 'charge']].groupby(
                ['accession', 'opt_global_cv_MS:1000889_peptidoform_sequence', 'charge']).size()
            counter.update(spectra_dict.to_dict())
        return counter

    def merge_psm_to_pep(self, mzTab_path, output_path, psm_chunksize=100000, pep_chunksize=100000):
        # load ms_runs
        ms_runs = self.extract_ms_runs(mzTab_path)

        spectra_count_dict = self.__get_spectra_count(
            mzTab_path, psm_chunksize)

        # header swith
        header = True
        if self.experiment_type != 'lfq':
            peps = self.skip_and_load_csv(mzTab_path, 'PEH', sep='\t', usecols=[
                                          "retention_time", 'accession', "opt_global_cv_MS:1000889_peptidoform_sequence", "best_search_engine_score[1]", "charge"], chunksize=pep_chunksize)
            for pep in peps:
                pep['retention_time'] = pep['retention_time'].round(1)
                pep.drop_duplicates(subset=['accession', 'retention_time', 'charge',
                                    'opt_global_cv_MS:1000889_peptidoform_sequence'], inplace=True)
                psms = self.skip_and_load_csv(
                    mzTab_path, 'PSH', sep='\t', usecols=self._tmt_usecols, chunksize=psm_chunksize)
                for psm in psms:
                    psm['scan_number'] = psm['spectra_ref'].apply(
                        lambda x: re.findall(r'scan=(\d+)', x)[0])
                    psm['spectra_ref'] = psm['spectra_ref'].apply(
                        lambda x: ms_runs[x.split(":")[0]])
                    psm['spectral_count'] = psm[['accession', 'opt_global_cv_MS:1000889_peptidoform_sequence', 'charge']].apply(
                        lambda row: spectra_count_dict[(row['accession'], row['opt_global_cv_MS:1000889_peptidoform_sequence'], row['charge'])], axis=1)
                    psm['retention_time'] = psm['retention_time'].round(1)
                    psm.drop_duplicates(subset=['accession', 'retention_time', 'charge',
                                        'opt_global_cv_MS:1000889_peptidoform_sequence'], inplace=True)

                    res = pep.merge(psm, on=['accession', 'retention_time', 'charge',
                                    'opt_global_cv_MS:1000889_peptidoform_sequence'], how='left')
                    res.dropna(subset=['spectra_ref'], inplace=True)
                    if header != True:
                        res.to_csv(output_path, mode='a+',
                                   index=False, header=False)
                    else:
                        header = False
                        res.to_csv(output_path, mode='a+', index=False)
        else:
            peps = self.skip_and_load_csv(mzTab_path, 'PEH', sep='\t', usecols=[
                                          'spectra_ref', 'accession', "opt_global_cv_MS:1000889_peptidoform_sequence", "best_search_engine_score[1]", "charge"], chunksize=pep_chunksize)
            for pep in peps:
                pep['spectra_ref'] = pep['spectra_ref'].apply(
                    lambda x: ms_runs[x.split(":")[0]])
                pep.drop_duplicates(subset=['accession', 'spectra_ref', 'charge',
                                    'opt_global_cv_MS:1000889_peptidoform_sequence'], inplace=True)
                psms = self.skip_and_load_csv(
                    mzTab_path, 'PSH', sep='\t', usecols=self._lfq_usecols, chunksize=psm_chunksize)
                for psm in psms:
                    psm['scan_number'] = psm['spectra_ref'].apply(
                        lambda x: re.findall(r'scan=(\d+)', x)[0])
                    psm['spectra_ref'] = psm['spectra_ref'].apply(
                        lambda x: ms_runs[x.split(":")[0]])
                    psm['spectral_count'] = psm[['accession', 'opt_global_cv_MS:1000889_peptidoform_sequence', 'charge']].apply(
                        lambda row: spectra_count_dict[(row['accession'], row['opt_global_cv_MS:1000889_peptidoform_sequence'], row['charge'])], axis=1)
                    psm['retention_time'] = psm['retention_time'].round(1)
                    psm.drop_duplicates(subset=['accession', 'spectra_ref', 'charge',
                                        'opt_global_cv_MS:1000889_peptidoform_sequence'], inplace=True)
                    res = pep.merge(psm, on=['accession', 'spectra_ref', 'charge',
                                    'opt_global_cv_MS:1000889_peptidoform_sequence'], how='left')
                    res.dropna(subset=['retention_time'], inplace=True)
                    if header != True:
                        res.to_csv(output_path, mode='a+',
                                   index=False, header=False)
                    else:
                        header = False
                        res.to_csv(output_path, mode='a+', index=False)

    def merge_pep_to_msstats_in(self, pep_path, msstats_path, output_path, msstats_chunksize=100000, pep_chunksize=100000):
        msstats_ins = pd.read_csv(msstats_path, chunksize=msstats_chunksize)
        header = True
        for msstats_in in msstats_ins:
            msstats_in['Reference'] = msstats_in['Reference'].apply(
                lambda x: x.split(".")[0])
            if 'FragmentIon' not in msstats_in.columns.to_list():
                msstats_in.loc[:, 'FragmentIon'] = 'Unknown'
            if 'IsotopeLabelType' not in msstats_in.columns.to_list():
                msstats_in.loc[:, 'IsotopeLabelType'] = 'L'

            if self.experiment_type != 'lfq':
                msstats_in["Channel"] = msstats_in["Channel"].apply(
                    lambda row: tmt[self.experiment_type][row-1])
                msstats_in["RetentionTime"] = msstats_in["RetentionTime"].round(
                    1)
                msstats_in = msstats_in[self._tmt_msstats_usecols]
                # msstats_in.drop_duplicates(subset=['ProteinName','RetentionTime','PeptideSequence','Charge'],inplace=True)
                peps = pd.read_csv(pep_path, chunksize=pep_chunksize)
                for pep in peps:
                    res = pd.merge(msstats_in, pep, left_on=['ProteinName', 'RetentionTime', 'PeptideSequence', 'Charge'], right_on=[
                                   'accession', 'retention_time', 'opt_global_cv_MS:1000889_peptidoform_sequence', 'charge'], how='left')
                    res.drop(['accession', 'retention_time',
                             'opt_global_cv_MS:1000889_peptidoform_sequence', 'charge'], axis=1, inplace=True)
                    res.dropna(subset=['spectra_ref'], inplace=True)
                    if header != True:
                        res.to_csv(output_path, mode='a+',
                                   index=False, header=False)
                    else:
                        header = False
                        res.to_csv(output_path, mode='a+', index=False)

            else:
                msstats_in = msstats_in[self._lfq_msstats_usecols]
                msstats_in.loc[:, 'Channel'] = 'label free sample'
                # msstats_in.drop_duplicates(subset=['ProteinName','Reference','PeptideSequence','PrecursorCharge'],inplace=True)
                peps = pd.read_csv(pep_path, chunksize=pep_chunksize)
                for pep in peps:
                    res = pd.merge(msstats_in, pep, left_on=['ProteinName', 'Reference', 'PeptideSequence', 'PrecursorCharge'], right_on=[
                                   'accession', 'spectra_ref', 'opt_global_cv_MS:1000889_peptidoform_sequence', 'charge'], how='left')
                    res.drop(['accession', 'spectra_ref',
                             'opt_global_cv_MS:1000889_peptidoform_sequence', 'charge'], axis=1, inplace=True)
                    res.dropna(subset=['retention_time'], inplace=True)
                    if header != True:
                        res.to_csv(output_path, mode='a+',
                                   index=False, header=False)
                    else:
                        header = False
                        res.to_csv(output_path, mode='a+', index=False)

    def merge_sdrf_to_msstats_in(self, sdrf_path, msstats_path, output_path, msstats_chunksize=100000):
        msstats_ins = pd.read_csv(msstats_path, chunksize=msstats_chunksize)
        sdrf = pd.read_csv(sdrf_path, sep='\t')
        factor = "".join(
            filter(lambda x: x.startswith('factor'), sdrf.columns))
        sdrf = sdrf[['comment[data file]', 'source name', factor,
                     'comment[fraction identifier]', 'comment[label]']]
        sdrf['comment[data file]'] = sdrf['comment[data file]'].apply(
            lambda x: x.split(".")[0])
        header = True
        if self.experiment_type != 'lfq':
            for msstats_in in msstats_ins:
                res = pd.merge(msstats_in, sdrf, left_on=['spectra_ref', "Channel"], right_on=[
                               'comment[data file]', 'comment[label]'], how='left')
                res.drop(['comment[data file]', 'comment[label]'],
                         axis=1, inplace=True)
                res.rename(columns=self._map_tmt, inplace=True)
                if header != True:
                    res.to_csv(output_path, mode='a+',
                               index=False, header=False)
                else:
                    header = False
                    res.to_csv(output_path, mode='a+', index=False)
        else:
            for msstats_in in msstats_ins:
                res = pd.merge(msstats_in, sdrf, left_on=['Reference'], right_on=[
                               'comment[data file]'], how='left')
                res.drop(['comment[data file]', 'comment[label]'],
                         axis=1, inplace=True)
                res.rename(columns=self._map_lfq, inplace=True)
                if header != True:
                    res.to_csv(output_path, mode='a+',
                               index=False, header=False)
                else:
                    header = False
                    res.to_csv(output_path, mode='a+', index=False)

    def skip_and_load_csv(self, fle, header, **kwargs):
        if self._mzTab_file == fle and self._psm_pos != None and header == 'PSH':
            return self.__load_second(fle, header, **kwargs)
        if self._mzTab_file == fle and self._pep_pos != None and header == 'PEH':
            return self.__load_second(fle, header, **kwargs)
        if self._mzTab_file == fle and self._prt_pos != None and header == 'PRH':
            return self.__load_second(fle, header, **kwargs)
        fle_len = self.__extract_len(fle, header)
        if os.stat(fle).st_size == 0:
            raise ValueError("File is empty")
        f = open(fle)
        pos = 0
        line = f.readline()
        while line.split("\t")[0] != header:
            pos = f.tell()
            line = f.readline()
        f.seek(pos)
        self.__set_table_config(header, fle_len, pos, fle)
        return pd.read_csv(f, nrows=fle_len, **kwargs)

    # extract ms runs
    def extract_ms_runs(self, fle):
        if os.stat(fle).st_size == 0:
            raise ValueError("File is empty")
        f = codecs.open(fle, 'r', 'utf-8')
        line = f.readline()
        ms_runs = {}
        while line.split("\t")[0] == 'MTD':
            if line.split("\t")[1].split("-")[-1] == 'location':
                ms_runs[line.split("\t")[1].split(
                    "-")[0]] = line.split("\t")[2].split("//")[-1].split(".")[0]
            line = f.readline()
        f.close()
        return ms_runs

    # optional

    def extract_optional_spectrum_featrue(self, res_df, mzml_directory):

        self.mzml_directory = mzml_directory
        res_df[['mz', 'array_intensity', 'num_peaks']] = res_df[['reference_file_name', 'scan_number']].apply(
            lambda x: self._map_spectrum_mz(x['reference_file_name'], x['scan_number']), axis=1, result_type="expand")

        return res_df

    def _map_spectrum_mz(self, mz_path, scan):
        scan = "controllerType=0 controllerNumber=1 scan=" + str(scan)
        mz_path = self.mzml_directory + '/' + mz_path + '.mzml'
        spectral = mzml.MzML(mz_path).get_by_id(scan)
        mz_array = spectral["m/z array"]
        array_intensity = spectral['intensity array']
        return mz_array, array_intensity, len(mz_array)

    def __split_start_or_end(self, value):
        if ',' in str(value):
            return list(map(int, value.split(',')))
        elif value is np.nan:
            return []
        else:
            return [int(value)]

    def convert_to_parquet(self, res_path):

        res = pd.read_csv(res_path)
        res['sequence'] = res['sequence'].fillna("").astype(str)
        res['protein_accessions'] = res['protein_accessions'].str.split('|')
        res['unique'] = res['unique'].fillna(0).astype(int)
        res['modifications'] = res['modifications'].fillna("").str.split(',')
        res['retention_time'] = res['retention_time'].fillna(0.0)
        res['charge'] = res['charge'].fillna(0).astype(int)
        res['exp_mass_to_charge'] = res['exp_mass_to_charge'].fillna(0.0)
        res['calc_mass_to_charge'] = res['calc_mass_to_charge'].fillna(0.0)
        res['posterior_error_probability'] = res['posterior_error_probability'].fillna(
            0.0)
        res['global_qvalue'] = res['global_qvalue'].fillna(0.0)
        res['is_decoy'] = res['is_decoy'].fillna(0).astype(int)
        res['protein_start_positions'] = res['protein_start_positions'].apply(
            self.__split_start_or_end)
        res['protein_end_positions'] = res['protein_end_positions'].apply(
            self.__split_start_or_end)

        res['best_id_score'] = res['best_id_score'].apply(
            lambda x: "\"OpenMS:Best_PSM_Score\"" + ": " + str(x))
        res['fraction'] = res['fraction'].fillna(1).astype(str)
        res['biological_replicate'] = res['biological_replicate'].fillna(
            1).astype(str)
        res['run'] = res['run'].fillna(1).astype(str)
        res['scan_number'] = res['scan_number'].astype(int).astype(str)

        # optional
        optional_list = ['gene_accessions', 'gene_names', 'id_scores']
        optional_float = ['mz', 'intensity_array']
        for op_l in optional_list:
            if op_l not in res.columns:
                res.loc[:, op_l] = "Unknow"
                res[op_l] = res[op_l].apply(lambda x: [x])

        for op_f in optional_float:
            if op_f not in res.columns:
                res.loc[:, op_f] = 0.0
                res[op_f] = res[op_f].apply(lambda x: [x])

        if 'consensus_support' not in res.columns:
            res.loc[:, 'consensus_support'] = 0.0
        else:
            res['consensus_support'] = res['consensus_support'].fillna(0.0)

        if 'num_peaks' not in res.columns:
            res.loc[:, 'num_peaks'] = 0

        return pa.Table.from_pandas(res, schema=self.schema)

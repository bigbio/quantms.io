import os
import pyarrow as pa
import pandas as pd
import numpy as np
import re
import codecs
from quantms_io.core.mztab import fetch_modifications_from_mztab_line
from quantms_io.utils.pride_utils import clean_peptidoform_sequence, get_petidoform_msstats_notation

'''
example
Convert = FeatureConvertor('lfq',schema)
Convert.merge_psm_to_pep("lfq2\PXD002854-serum.sdrf_openms_design_openms.mzTab",'res1.txt')
Convert.merge_pep_and_sdrf_to_msstats_in("res1.txt","lfq2\PXD002854-serum.sdrf_openms_design_msstats_in.csv","lfq2\PXD002854-serum.sdrf.tsv","result_lfq.csv")
parquet = Convert.convert_to_parquet("result_lfq.csv")
'''
from quantms_io.utils.constants import TMT_CHANNELS, ITRAQ_CHANNEL

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
        # load psms columns
        self._psms_columns = None
        #load pep columns
        self._pep_columns = None
        #ms_runs
        self._ms_runs = None
        # msstats_in usecols
        self._tmt_msstats_usecols = ['ProteinName', 'PeptideSequence', 'Channel', 'Run', 'BioReplicate',
                                     'Intensity', 'FragmentIon', 'IsotopeLabelType', 'Charge', 'Reference','protein_global_qvalue','sequence']
        self._lfq_msstats_usecols = ['ProteinName', 'PeptideSequence', 'Channel', 'Reference',
                                     'Run', 'BioReplicate', 'Intensity', 'FragmentIon', 'IsotopeLabelType', 'PrecursorCharge','protein_global_qvalue','sequence']
        # score names
        self._score_names = None  # dict
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
            'source name': 'sample_accession',
            'comment[fraction identifier]': 'fraction',
            'factor value[organism part]': 'condition',
        }
        self._map_tmt = {
            'ProteinName': 'protein_accessions',
            'PeptideSequence': 'peptidoform',
            'Reference':  'reference_file_name',
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
            'factor value[organism part]': 'condition',
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
    
    # extract pep columns
    def __extract_pep_columns(self,file):
        if os.stat(file).st_size == 0:
            raise ValueError("File is empty")
        f = open(file)
        line = f.readline()
        while not line.startswith("PEH"):
            line = f.readline()
        self._pep_columns = line.split('\n')[0].split('\t')

    # extract csv len
    def __extract_len(self, fle, header):
        map_tag = {
            "PSH": 'PSM',
            "PEH": 'PEP',
            'PRH': 'PRT'
        }
        if os.stat(fle).st_size == 0:
            raise ValueError("File is empty")
        f = open(fle)
        pos = 0
        line = f.readline()
        while line.split("\t")[0] != header:
            pos = f.tell()
            line = f.readline()
            
        if header == 'PSH':
            self._psms_columns = line.split('\n')[0].split('\t')
        
        line = f.readline()
        fle_len = 0
        while line.split("\t")[0] == map_tag[header]:
            fle_len += 1
            line = f.readline()
        f.close()
        return fle_len,pos
    
    def skip_and_load_csv(self, fle, header, **kwargs):
        '''
        file: mzTab file
        header: table tag #PSH PRH PEH
        '''
        if self._mzTab_file == fle and self._psm_pos != None and header == 'PSH':
            return self.__load_second(fle, header, **kwargs)
        if self._mzTab_file == fle and self._pep_pos != None and header == 'PEH':
            return self.__load_second(fle, header, **kwargs)
        if self._mzTab_file == fle and self._prt_pos != None and header == 'PRH':
            return self.__load_second(fle, header, **kwargs)
        fle_len,pos = self.__extract_len(fle, header)
        if os.stat(fle).st_size == 0:
            raise ValueError("File is empty")
        f = open(fle)
        f.seek(pos)
        self.__set_table_config(header,fle_len, pos, fle)
        return pd.read_csv(f, nrows=fle_len, **kwargs)

    def __get_spectra_count(self, mzTab_path,psm_chunksize):
        '''
        mzTab_path: mzTab file path
        psm_chunksize: the large of chunk
        return: a dict about piptie numbers
        '''
        from collections import Counter
        counter = Counter()
        psms = self.skip_and_load_csv(
            mzTab_path, 'PSH', sep='\t', chunksize=psm_chunksize)
        for psm in psms:
            psm['spectra_ref'] = psm['spectra_ref'].apply(lambda x: self._ms_runs[x.split(":")[0]])
            spectra_dict = psm[['opt_global_cv_MS:1000889_peptidoform_sequence', 'charge','spectra_ref']].groupby(
                ['opt_global_cv_MS:1000889_peptidoform_sequence', 'charge','spectra_ref']).size()
            counter.update(spectra_dict.to_dict())
        return counter

    def __get_protein_map(self, mzTab_path):
        '''
        return: a dict about protein score
        '''
        prt = self.skip_and_load_csv(mzTab_path, 'PRH', sep='\t', usecols=[
                                     'accession', 'best_search_engine_score[1]'])
        prt_score = prt.groupby('accession').max()
        protein_map = prt_score.to_dict()['best_search_engine_score[1]']
        return protein_map

    def __get_modifications(self, fle):
        '''
        return: a dict about modifications
        '''
        if os.stat(fle).st_size == 0:
            raise ValueError("File is empty")
        f = codecs.open(fle, 'r', 'utf-8')
        line = f.readline()
        mod_dict = {}
        while line.split("\t")[0] == 'MTD':
            if "_mod[" in line:
                mod_dict = fetch_modifications_from_mztab_line(line, mod_dict)
            line = f.readline()
        f.close()
        return mod_dict

    def __get_score_names(self, fle):
        '''
        return: a dict about search engine
        '''
        if os.stat(fle).st_size == 0:
            raise ValueError("File is empty")
        f = codecs.open(fle, 'r', 'utf-8')
        line = f.readline()
        score_names = {}
        while line.split("\t")[0] == 'MTD':
            if 'search_engine_score' in line:
                score_values = line.replace(
                    "[", "").replace("]", "").split(",")
                score_name = score_values[2].strip()
                if ":" in score_name:
                    score_name = "'{}'".format(score_name)
                if 'peptide' in line.split('\t')[1]:
                    score_names['peptide_score'] = score_name
                elif 'protein' in line.split('\t')[1]:
                    score_names['protein_score'] = score_name
                else:
                    score_names['psm_score'] = score_name

            line = f.readline()
        f.close()
        return score_names

    def __handle_protein_map(self, protein_map, key):
        '''
        map protein score from accession
        '''
        if key not in protein_map.keys():
            keys = key.split(',')
            for k in keys:
                if k in protein_map.keys():
                    return protein_map[k]
            return None
        else:
            return protein_map[key]
        
    def _extract_from_pep(self,mzTab_path):
        '''
        return: dict about pep_msg
        '''
        self.__extract_pep_columns(mzTab_path)
        pep_usecols = ['opt_global_cv_MS:1000889_peptidoform_sequence','charge','best_search_engine_score[1]','spectra_ref']
        live_cols = [col for col in pep_usecols if col in self._pep_columns]
        not_cols = [col for col in pep_usecols if col not in live_cols]
        if 'opt_global_cv_MS:1000889_peptidoform_sequence' in not_cols:
            if 'sequence' in self._pep_columns and 'modifications' in self._pep_columns:
                live_cols.append('sequence')
                live_cols.append('modifications')
            else:
                raise Exception("The peptide table don't have opt_global_cv_MS:1000889_peptidoform_sequence columns")
        if 'charge' in not_cols or 'best_search_engine_score[1]' in not_cols:
                raise Exception("The peptide table don't have best_search_engine_score[1] or charge columns")
                
        pep = self.skip_and_load_csv(mzTab_path,'PEH',sep='\t',usecols=live_cols)
        
        # check opt_global_cv_MS:1000889_peptidoform_sequence
        if 'opt_global_cv_MS:1000889_peptidoform_sequence' not in pep.columns:
            modifications = self.__get_modifications(mzTab_path)
            pep['opt_global_cv_MS:1000889_peptidoform_sequence'] = pep[['modifications', 'sequence']].apply(lambda row: get_petidoform_msstats_notation(row['sequence'], row['modifications'], modifications), axis=1)
        
        # check spectra_ref
        if 'spectra_ref' not in pep.columns:
            pep.loc[:,'scan_number'] = None
            pep.loc[:,'spectra_ref'] = None
        else:
            pep['scan_number'] = pep['spectra_ref'].apply(lambda x: re.findall(r'scan=(\d+)',x)[0])
            pep['spectra_ref'] = pep['spectra_ref'].apply(lambda x:self._ms_runs[x.split(":")[0]])
        pep_msg = pep.groupby(['opt_global_cv_MS:1000889_peptidoform_sequence','charge']).max()
        pep_msg['pep_msg'] = pep_msg[['best_search_engine_score[1]','spectra_ref','scan_number']].apply(lambda row: [row['best_search_engine_score[1]'],row['spectra_ref'],row['scan_number']],axis=1)
        
        map_dict = pep_msg.to_dict()['pep_msg']
        return map_dict
    
    def _extract_from_psm_to_pep_msg(self,mzTab_path,map_dict):
        '''
        return dict about pep and psm msg
        '''
        psm = self.skip_and_load_csv(mzTab_path,'PSH',sep='\t')
        modifications = self.__get_modifications(mzTab_path)
        #psm['spectra_ref'] = psm['spectra_ref'].apply(lambda x:ms_runs[x.split(":")[0]])
        if 'opt_global_cv_MS:1000889_peptidoform_sequence' not in psm.columns:
            psm['opt_global_cv_MS:1000889_peptidoform_sequence'] = psm[['modifications', 'sequence']].apply(lambda row: get_petidoform_msstats_notation(row['sequence'], row['modifications'], modifications), axis=1)
        for key,df in psm.groupby(['opt_global_cv_MS:1000889_peptidoform_sequence','charge']):
            if key not in map_dict.keys():
                map_dict[key] = [None,None,None]
            
            if pd.isna(map_dict[key][1]):
                df['scan_number'] = df['spectra_ref'].apply(lambda x: re.findall(r'scan=(\d+)',x)[0])
                df['spectra_ref'] = df['spectra_ref'].apply(lambda x:self._ms_runs[x.split(":")[0]])
                if 'opt_global_q-value_score' in df.columns:
                    map_dict[key][1] = df[df['opt_global_q-value_score'] == df['opt_global_q-value_score'].max()]['spectra_ref'].values[0]
                    map_dict[key][2] = df[df['opt_global_q-value_score'] == df['opt_global_q-value_score'].max()]['scan_number'].values[0]
                elif 'search_engine_score[1]' in df.columns:
                    map_dict[key][1] = df[df['search_engine_score[1]'] == df['search_engine_score[1]'].max()]['spectra_ref'].values[0]
                    map_dict[key][2] = df[df['search_engine_score[1]'] == df['search_engine_score[1]'].max()]['scan_number'].values[0]
                else:
                    raise Exception("The psm table don't have opt_global_q-value_score or search_engine_score[1] columns")
            else:
                df['spectra_ref'] = df['spectra_ref'].apply(lambda x:self._ms_runs[x.split(":")[0]])
            map_dict[key].append(df['start'].values[0])
            map_dict[key].append(df['end'].values[0])
            map_dict[key].append(df['unique'].values[0])
            map_dict[key].append(df['modifications'].values[0])
            if 'opt_global_Posterior_Error_Probability_score' in df.columns:
                map_dict[key].append(df['opt_global_Posterior_Error_Probability_score'].values[0])
            else:
                map_dict[key].append(None)
            if 'opt_global_cv_MS:1002217_decoy_peptide' in df.columns:
                map_dict[key].append(df['opt_global_cv_MS:1002217_decoy_peptide'].values[0])
            else:
                map_dict[key].append(None)
            if map_dict[key][1] not in df['spectra_ref'].values:
                map_dict[key].append(df["calc_mass_to_charge"].values[0])
                map_dict[key].append(None)
            else:
                map_dict[key].append(df[df['spectra_ref']==map_dict[key][1]]["calc_mass_to_charge"].values[0])
                map_dict[key].append(df[df['spectra_ref']==map_dict[key][1]]["exp_mass_to_charge"].values[0])
                
        return map_dict
        
    def _extract_psm_pep_msg(self, mzTab_path):
        '''
        mzTab_path: mzTab file path
        output_path: output file path (csv)
        psm_chunksize: the large of psm chunk
        pep_chunksize: the large of pep chunk
        '''
        # load ms_runs
        self._ms_runs = self.extract_ms_runs(mzTab_path)
        self._score_names = self.__get_score_names(mzTab_path)
        map_dict = self._extract_from_pep(mzTab_path)
        map_dict = self._extract_from_psm_to_pep_msg(mzTab_path,map_dict)
        
        return map_dict
    
    def _map_msstats_in(self,msstats_in,map_dict,spectra_count_dict):
        '''
        map key: PeptideSequence,charge
        '''
        if self.experiment_type == 'LFQ':
            msstats_in["spectral_count"] = msstats_in[['PeptideSequence','PrecursorCharge','Reference']].apply(lambda row: spectra_count_dict[(row['PeptideSequence'], row['PrecursorCharge'],row['Reference'])], axis=1)
        else:
            msstats_in["spectral_count"] = msstats_in[['PeptideSequence','Charge','Reference']].apply(lambda row: spectra_count_dict[(row['PeptideSequence'], row['Charge'],row['Reference'])], axis=1)
        map_features = ["global_qvalue","best_psm_reference_file_name","best_psm_scan_number","protein_start_positions","protein_end_positions",'unique','modifications',"posterior_error_probability","is_decoy","calc_mass_to_charge","exp_mass_to_charge"]
        for i,feature in enumerate(map_features):
            if self.experiment_type == 'LFQ':
                msstats_in.loc[:,feature] = msstats_in[['PeptideSequence','PrecursorCharge']].apply(lambda row:map_dict[(row['PeptideSequence'],row['PrecursorCharge'])][i],axis=1)
            else:
                msstats_in.loc[:,feature] = msstats_in[['PeptideSequence','Charge']].apply(lambda row:map_dict[(row['PeptideSequence'],row['Charge'])][i],axis=1)   
        peptide_score_name = self._score_names['peptide_score']
        msstats_in["id_scores"] = msstats_in[["global_qvalue","posterior_error_probability"]].apply(lambda row: peptide_score_name + ':' + str(row['global_qvalue']) +',' + 'Best PSM PEP:' + str(row['posterior_error_probability']),axis=1)
        return msstats_in
        
    def merge_mzTab_and_sdrf_to_msstats_in(self, mzTab_path, msstats_path, sdrf_path, output_path, msstats_chunksize=500000):
        '''
        pep_path: the output file of function merge_psm_to_pep 
        msstats_path: msstats_in file path
        sdrf_path: sdrf file path
        output_path: output file path(csv)
        msstats_chunksize: the large of msstats chunk
        pep_chunksize: the large of pep chunk
        '''
        protein_map = self.__get_protein_map(mzTab_path)
        map_dict = self._extract_psm_pep_msg(mzTab_path)
        msstats_ins = pd.read_csv(msstats_path, chunksize=msstats_chunksize)
        spectra_count_dict = self.__get_spectra_count(mzTab_path,500000)
        header = True
        for msstats_in in msstats_ins:
            msstats_in['Reference'] = msstats_in['Reference'].apply(
                lambda x: x.split(".")[0])
            msstats_in['protein_global_qvalue'] =  msstats_in['ProteinName'].apply(lambda x: self.__handle_protein_map(protein_map, x))
            msstats_in['sequence'] = msstats_in['PeptideSequence'].apply(lambda x: clean_peptidoform_sequence(x))
            if self.experiment_type != 'LFQ':
                no_tmt_usecols = [
                    col for col in self._tmt_msstats_usecols if col not in msstats_in.columns]
                for col in no_tmt_usecols:
                    if 'IsotopeLabelType' == col:
                        msstats_in.loc[:, col] = 'L'
                    else:
                        msstats_in.loc[:, col] = None
                if 'TMT' in self.experiment_type:
                    msstats_in["Channel"] = msstats_in["Channel"].apply(
                        lambda row: TMT_CHANNELS[self.experiment_type][row-1])
                else:
                    msstats_in["Channel"] = msstats_in["Channel"].apply(
                        lambda row: ITRAQ_CHANNEL[self.experiment_type][row-1])
                if 'RetentionTime' in msstats_in.columns:
                    self._tmt_msstats_usecols.append('RetentionTime')
                msstats_in = msstats_in[self._tmt_msstats_usecols]
                self._tmt_msstats_usecols.remove('RetentionTime')
                msstats_in = self._map_msstats_in(msstats_in,map_dict,spectra_count_dict)
                self._merge_sdrf_to_msstats_in(sdrf_path, msstats_in, output_path, header, msstats_chunksize)
                header = False
            else:
                no_lfq_usecols = [
                    col for col in self._lfq_msstats_usecols if col not in msstats_in.columns]
                for col in no_lfq_usecols:
                    if col == 'Channel':
                        msstats_in.loc[:, col] = 'label free sample'
                    else:
                        msstats_in.loc[:, col] = None
                msstats_in = msstats_in[self._lfq_msstats_usecols]
                msstats_in = self._map_msstats_in(msstats_in,map_dict,spectra_count_dict)
                self._merge_sdrf_to_msstats_in(sdrf_path,msstats_in, output_path, header, msstats_chunksize)
                header = False

    def _merge_sdrf_to_msstats_in(self, sdrf_path, msstats_in, output_path, header, msstats_chunksize):
        '''
        sdrf_path:
        msstats_in: dataframe of msstats
        '''
        sdrf = pd.read_csv(sdrf_path, sep='\t')
        factor = "".join(
            filter(lambda x: x.startswith('factor'), sdrf.columns))
        sdrf = sdrf[['comment[data file]', 'source name', factor,
                     'comment[fraction identifier]', 'comment[label]']]
        sdrf['comment[data file]'] = sdrf['comment[data file]'].apply(
            lambda x: x.split(".")[0])
        if self.experiment_type != 'LFQ':
            res = pd.merge(msstats_in, sdrf, left_on=['Reference', "Channel"], right_on=[
                           'comment[data file]', 'comment[label]'], how='left')
            res.drop(['comment[data file]', 'comment[label]'],
                     axis=1, inplace=True)
            res.rename(columns=self._map_tmt, inplace=True)
            if header != True:
                res.to_csv(output_path, mode='a+', index=False, header=False,sep='\t')
            else:
                header = False
                res.to_csv(output_path, mode='a+', index=False, sep='\t')
        else:
            res = pd.merge(msstats_in, sdrf, left_on=['Reference'], right_on=[
                           'comment[data file]'], how='left')
            res.drop(['comment[data file]', 'comment[label]'],
                     axis=1, inplace=True)
            res.rename(columns=self._map_lfq, inplace=True)
            if header != True:
                res.to_csv(output_path, mode='a+', index=False, header=False,sep='\t')
            else:
                header = False
                res.to_csv(output_path, mode='a+', index=False,sep='\t')

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
        # TODO: Replace by OpenMSHandler function
        # scan = "controllerType=0 controllerNumber=1 scan=" + str(scan)
        # mz_path = self.mzml_directory + '/' + mz_path + '.mzml'
        # spectral = mzml.MzML(mz_path).get_by_id(scan)
        # mz_array = spectral["m/z array"]
        # array_intensity = spectral['intensity array']
        return [], [], 0

    def __split_start_or_end(self, value):
        if pd.isna(value):
            return pd.NA
        elif ',' in str(value):
            return list(map(int, value.split(',')))
        elif value is np.nan:
            return None
        else:
            return [int(value)]

    def convert_to_parquet(self, res_path):
        '''
        return: parquet
        '''
        res = pd.read_csv(res_path,sep='\t')
        res['id_scores'] = res['id_scores'].str.split(',')
        res['sequence'] = res['sequence'].astype(str)
        res['protein_accessions'] = res['protein_accessions'].str.split(",")
        res['protein_start_positions'] = res['protein_start_positions'].apply(
            self.__split_start_or_end).to_list()
        res['protein_end_positions'] = res['protein_end_positions'].apply(
            self.__split_start_or_end).to_list()
        res['protein_global_qvalue'] = res['protein_global_qvalue'].astype(float)
        res['unique'] = res['unique'].map(lambda x: pd.NA if pd.isna(x) else int(x)).astype('Int32')
        res['modifications'] = res['modifications'].str.split(",")
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
        res['best_psm_scan_number'] = res['best_psm_scan_number'].astype(int).astype(str)
        
        if 'retention_time' in res.columns:
            res['retention_time'] = res['retention_time'].astype(float)
        else:
            res.loc[:,'retention_time'] = None
        
        res.loc[:,'gene_accessions'] = None
        res.loc[:,'gene_names'] = None
        res.loc[:,'num_peaks'] = None
        res.loc[:,'mz'] = None
        res.loc[:,'intensity_array'] = None
        return pa.Table.from_pandas(res, schema=self.schema)

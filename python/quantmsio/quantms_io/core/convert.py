import os
import pyarrow as pa
import pandas as pd
import numpy as np
from pyteomics import mzml
import codecs
from quantms_io.core.mztab import fetch_modifications_from_mztab_line, get_petidoform_msstats_notation
'''
example
Convert = FeatureConvertor('lfq',schema)
Convert.merge_psm_to_pep("lfq2\PXD002854-serum.sdrf_openms_design_openms.mzTab",'res1.txt')
Convert.merge_pep_and_sdrf_to_msstats_in("res1.txt","lfq2\PXD002854-serum.sdrf_openms_design_msstats_in.csv","lfq2\PXD002854-serum.sdrf.tsv","result_lfq.csv")
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
        # load psms columns
        self._psms_columns = None
        # ms_runs
        self._ms_runs = None
        # psm usecols
        self._psm_usecols = ['sequence', 'accession', 'unique', 'modifications', 'retention_time', 'charge', 'exp_mass_to_charge', 'calc_mass_to_charge', 'start', 'end',
                             'opt_global_Posterior_Error_Probability_score', 'opt_global_cv_MS:1000889_peptidoform_sequence', 'opt_global_cv_MS:1002217_decoy_peptide', 'opt_global_consensus_support', 'spectra_ref']
        # msstats_in usecols
        self._tmt_msstats_usecols = ['ProteinName', 'PeptideSequence', 'Channel', 'Run', 'BioReplicate',
                                     'Intensity', 'FragmentIon', 'IsotopeLabelType', 'RetentionTime', 'Charge', 'Reference']
        self._lfq_msstats_usecols = ['ProteinName', 'PeptideSequence', 'Channel', 'Reference',
                                     'Run', 'BioReplicate', 'Intensity', 'FragmentIon', 'IsotopeLabelType', 'PrecursorCharge']
        # score names
        self._score_names = None  # dict
        # modifications
        # self._modifications = None #dict
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
            'best_search_engine_score[1]': 'global_qvalue',
            'start': 'protein_start_positions',
            'end': 'protein_end_positions',
            'opt_global_Posterior_Error_Probability_score': 'posterior_error_probability',
            'opt_global_cv_MS:1002217_decoy_peptide': 'is_decoy',
            'opt_global_consensus_support': 'consensus_support',
            'source name': 'sample_accession',
            'comment[fraction identifier]': 'fraction',
            'factor value[organism part]': 'condition'
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
            'best_search_engine_score[1]': 'global_qvalue',
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
        if header == 'PSH':
            self._psms_columns = line.split('\n')[0].split('\t')
        line = f.readline()
        fle_len = 0
        while line.split("\t")[0] == map_tag[header]:
            fle_len += 1
            line = f.readline()
        f.close()
        return fle_len

    def __get_spectra_count(self, mzTab_path, psm_chunksize):
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
            psm['spectra_ref'] = psm['spectra_ref'].apply(
                lambda x: self._ms_runs[x.split(":")[0]])
            spectra_dict = psm[['opt_global_cv_MS:1000889_peptidoform_sequence', 'charge', 'spectra_ref']].groupby(
                ['opt_global_cv_MS:1000889_peptidoform_sequence', 'charge', 'spectra_ref']).size()
            counter.update(spectra_dict.to_dict())
        return counter

    def __get_protein_map(self, mzTab_path):
        '''
        return: a dict about protein score
        '''
        prt = self.skip_and_load_csv(mzTab_path, 'PRH', sep='\t', usecols=[
                                     'accession', 'best_search_engine_score[1]'])
        prt_score = prt.groupby('accession').mean()
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

    def merge_psm_to_pep(self, mzTab_path, output_path, psm_chunksize=100000, pep_chunksize=100000):
        '''
        mzTab_path: mzTab file path
        output_path: output file path (csv)
        psm_chunksize: the large of psm chunk
        pep_chunksize: the large of pep chunk
        '''
        # load ms_runs
        self._ms_runs = self.extract_ms_runs(mzTab_path)
        # load protein_map
        protein_map = self.__get_protein_map(mzTab_path)
        modifications = self.__get_modifications(mzTab_path)
        self._score_names = self.__get_score_names(mzTab_path)

        spectra_count_dict = self.__get_spectra_count(
            mzTab_path, psm_chunksize)
        be_usecols = [
            col for col in self._psm_usecols if col in self._psms_columns]
        no_usecols = [
            col for col in self._psm_usecols if col not in be_usecols]

        # header swith
        header = True
        peps = self.skip_and_load_csv(mzTab_path, 'PEH', sep='\t', usecols=[
                                      'spectra_ref', "opt_global_cv_MS:1000889_peptidoform_sequence", "best_search_engine_score[1]", "charge"], chunksize=pep_chunksize)
        for pep in peps:
            pep['spectra_ref'] = pep['spectra_ref'].apply(
                lambda x: self._ms_runs[x.split(":")[0]])
            psms = self.skip_and_load_csv(
                mzTab_path, 'PSH', sep='\t', usecols=be_usecols, chunksize=psm_chunksize)
            for psm in psms:
                for col in no_usecols:
                    if col == 'opt_global_cv_MS:1000889_peptidoform_sequence':
                        psm['opt_global_cv_MS:1000889_peptidoform_sequence'] = psm[['modifications', 'sequence']].apply(
                            lambda row: get_petidoform_msstats_notation(row['sequence'], row['modifications'], modifications), axis=1)
                    else:
                        psm.loc[:, col] = None
                psm['scan_number'] = psm['spectra_ref'].apply(
                    lambda x: re.findall(r'scan=(\d+)', x)[0])
                psm['spectra_ref'] = psm['spectra_ref'].apply(
                    lambda x: self._ms_runs[x.split(":")[0]])
                psm['protein_global_qvalue'] = psm['accession'].apply(
                    lambda x: self.__handle_protein_map(protein_map, x))
                psm['spectral_count'] = psm[['opt_global_cv_MS:1000889_peptidoform_sequence', 'charge', 'spectra_ref']].apply(
                    lambda row: spectra_count_dict[(row['opt_global_cv_MS:1000889_peptidoform_sequence'], row['charge'], row['spectra_ref'])], axis=1)
                psm.drop_duplicates(subset=[
                                    'spectra_ref', 'charge', 'opt_global_cv_MS:1000889_peptidoform_sequence'], inplace=True)
                res = pep.merge(psm, on=[
                                'spectra_ref', 'charge', 'opt_global_cv_MS:1000889_peptidoform_sequence'], how='left')
                res.dropna(subset=['retention_time'], inplace=True)
                if header != True:
                    res.to_csv(output_path, mode='a+',
                               index=False, header=False)
                else:
                    header = False
                    res.to_csv(output_path, mode='a+', index=False)

    def merge_pep_and_sdrf_to_msstats_in(self, pep_path, msstats_path, sdrf_path, output_path, msstats_chunksize=100000, pep_chunksize=100000):
        '''
        pep_path: the output file of function merge_psm_to_pep 
        msstats_path: msstats_in file path
        sdrf_path: sdrf file path
        output_path: output file path(csv)
        msstats_chunksize: the large of msstats chunk
        pep_chunksize: the large of pep chunk
        '''
        msstats_ins = pd.read_csv(msstats_path, chunksize=msstats_chunksize)
        header = True
        for msstats_in in msstats_ins:
            msstats_in['Reference'] = msstats_in['Reference'].apply(
                lambda x: x.split(".")[0])
            if self.experiment_type != 'lfq':
                no_tmt_usecols = [
                    col for col in self._tmt_msstats_usecols if col not in msstats_in.columns]
                for col in no_tmt_usecols:
                    msstats_in.loc[:, col] = None
                msstats_in["Channel"] = msstats_in["Channel"].apply(
                    lambda row: tmt[self.experiment_type][row-1])
                # msstats_in["RetentionTime"] = msstats_in["RetentionTime"].round(1)
                msstats_in = msstats_in[self._tmt_msstats_usecols]
                # msstats_in.drop_duplicates(subset=['ProteinName','RetentionTime','PeptideSequence','Charge'],inplace=True)
                peps = pd.read_csv(pep_path, chunksize=pep_chunksize)
                for pep in peps:
                    pep.drop_duplicates(subset=[
                                        'spectra_ref', 'charge', 'opt_global_cv_MS:1000889_peptidoform_sequence'], inplace=True)
                    res = pd.merge(msstats_in, pep, left_on=['Reference', 'PeptideSequence', 'Charge'], right_on=[
                                   'spectra_ref', 'opt_global_cv_MS:1000889_peptidoform_sequence', 'charge'], how='left')
                    res.drop(['accession', 'spectra_ref', 'opt_global_cv_MS:1000889_peptidoform_sequence',
                             'charge', 'retention_time'], axis=1, inplace=True)
                    res.dropna(subset=['start'], inplace=True)
                    self._merge_sdrf_to_msstats_in(
                        sdrf_path, res, output_path, header, msstats_chunksize)
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
                # msstats_in.drop_duplicates(subset=['ProteinName','Reference','PeptideSequence','PrecursorCharge'],inplace=True)
                peps = pd.read_csv(pep_path, chunksize=pep_chunksize)
                for pep in peps:
                    res = pd.merge(msstats_in, pep, left_on=['Reference', 'PeptideSequence', 'PrecursorCharge'], right_on=[
                                   'spectra_ref', 'opt_global_cv_MS:1000889_peptidoform_sequence', 'charge'], how='left')
                    res.drop(['accession', 'spectra_ref',
                             'opt_global_cv_MS:1000889_peptidoform_sequence', 'charge'], axis=1, inplace=True)
                    res.dropna(subset=['retention_time'], inplace=True)
                    self._merge_sdrf_to_msstats_in(
                        sdrf_path, res, output_path, header, msstats_chunksize)
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
        if self.experiment_type != 'lfq':
            res = pd.merge(msstats_in, sdrf, left_on=['Reference', "Channel"], right_on=[
                           'comment[data file]', 'comment[label]'], how='left')
            res.drop(['comment[data file]', 'comment[label]'],
                     axis=1, inplace=True)
            res.rename(columns=self._map_tmt, inplace=True)
            if header != True:
                res.to_csv(output_path, mode='a+', index=False, header=False)
            else:
                header = False
                res.to_csv(output_path, mode='a+', index=False)
        else:
            res = pd.merge(msstats_in, sdrf, left_on=['Reference'], right_on=[
                           'comment[data file]'], how='left')
            res.drop(['comment[data file]', 'comment[label]'],
                     axis=1, inplace=True)
            res.rename(columns=self._map_lfq, inplace=True)
            if header != True:
                res.to_csv(output_path, mode='a+', index=False, header=False)
            else:
                header = False
                res.to_csv(output_path, mode='a+', index=False)

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

    def query_accession(self, gene_name):
        Entrez.email = self.Email
        Entrez.tool = self.tool
        if gene_name == 'unknow':
            return ['unknow']
        q_name = gene_name + " AND Homo sapiens[porgn:__txid9606]"
        handle = Entrez.esearch(db="nucleotide", term=q_name)
        record = Entrez.read(handle)
        id_list = record["IdList"]
        accessions = get_accessions(id_list)
        return accessions

    def get_accessions(self, id_list):
        accessions = []
        for id_q in id_list:
            handle = Entrez.efetch(
                db="nucleotide", id=id_q, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            accessions.append(record.id)
        return accessions

    def extract_optional_gene_featrue(self, res_df, mzTab_path, gene_accession=False, Email=None, Tool=None):
        PRT = self.__skip_and_load_csv(mzTab_path, 'PRH', usecols=[
                                       'accession', 'description'])
        PRT['accession'] = PRT['accession'].apply(lambda x: x.split("|")[-1])
        PRT.loc[:, 'gene_names'] = PRT['description'].fillna('unkonw').apply(
            lambda x: re.findall(r'GN=(\w+)', x)[0]if re.findall(r'GN=(\w+)', x) else 'unknow')
        gene_df = protein[["accession", "gene_names"]].set_index(
            "accession").to_dict(orient='dict')["gene_names"]
        gene_df = pd.DataFrame(
            index=res.keys(), data=res.values()).reset_index()
        gene_df.columns = ['accession', 'gene_names']
        if gene_accession:
            self.email = Email
            self.tool = Tool
            gene_df.loc[:, 'gene_accession'] = gene_df['gene_names'].apply(
                lambda x: self.query_accession(x))
        res_df['accession'] = res_df['protein_accessions'].apply(
            lambda x: x[-1])
        res_df = res_df.merge(gene_df, on='accession', how='left')
        res_df.drop("accession", axis=1, inplace=True)
        return res_df

    def __split_start_or_end(self, value):
        if ',' in str(value):
            return list(map(int, value.split(',')))
        elif value is np.nan:
            return []
        else:
            return [int(value)]

    def convert_to_parquet(self, res_path):
        '''
        return: parquet
        '''
        res = pd.read_csv(res_path)
        res['sequence'] = res['sequence'].fillna("").astype(str)
        res['protein_accessions'] = res['protein_accessions'].str.split(",")
        res['protein_start_positions'] = res['protein_start_positions'].apply(
            self.__split_start_or_end).to_list()
        res['protein_end_positions'] = res['protein_end_positions'].apply(
            self.__split_start_or_end).to_list()
        res['protein_global_qvalue'] = res['protein_global_qvalue'].fillna(
            0.0).astype(float)
        res.loc[:, 'protein_best_id_score'] = res['protein_global_qvalue'].apply(
            lambda x: self._score_names['protein_score'] + ":" + str(round(float(x), 2)))
        res['unique'] = res['unique'].fillna(0)
        res['modifications'] = res['modifications'].fillna("").str.split(",")
        res['retention_time'] = res['retention_time'].fillna(0.0).astype(float)
        res['charge'] = res['charge'].fillna(0).astype(int)
        res['exp_mass_to_charge'] = res['exp_mass_to_charge'].fillna(
            0.0).astype(float)
        res['calc_mass_to_charge'] = res['calc_mass_to_charge'].fillna(
            0.0).astype(float)
        res['peptidoform'] = res['peptidoform'].fillna("")
        res['posterior_error_probability'] = res['posterior_error_probability'].fillna(
            0.0).astype(float)
        res['global_qvalue'] = res['global_qvalue'].fillna(0.0).astype(float)
        res['is_decoy'] = res['is_decoy'].fillna(0).astype(int)
        res.loc[:, 'best_id_score'] = res['global_qvalue'].apply(
            lambda x: self._score_names['peptide_score'] + ":" + str(round(float(x), 2)))
        res['intensity'] = res['intensity'].fillna(0.0).astype(float)
        res['spectral_count'] = res['spectral_count'].fillna(1).astype(int)
        res['sample_accession'] = res['sample_accession'].fillna("")
        res['condition'] = res['condition'].fillna("")
        res['fraction'] = res['fraction'].fillna("").astype(str)
        res['biological_replicate'] = res['biological_replicate'].fillna(
            "").astype(str)
        res['fragment_ion'] = res['fragment_ion'].fillna("").astype(str)
        res['isotope_label_type'] = res['isotope_label_type'].fillna("")
        res['run'] = res['run'].fillna("").astype(str)
        res['channel'] = res['channel'].fillna("")
        res['scan_number'] = res['scan_number'].astype(int).astype(str)
        res['consensus_support'] = res['consensus_support'].fillna(
            0.0).astype(float)

        # optional
        optional_list = ['gene_accessions', 'gene_names', 'id_scores']
        optional_float = ['mz', 'intensity_array']
        for op_l in optional_list:
            if op_l not in res.columns:
                res.loc[:, op_l] = ""
                res[op_l] = res[op_l].apply(lambda x: [x])

        for op_f in optional_float:
            if op_f not in res.columns:
                res.loc[:, op_f] = 0.0
                res[op_f] = res[op_f].apply(lambda x: [x])

        if 'num_peaks' not in res.columns:
            res.loc[:, 'num_peaks'] = 0

        return pa.Table.from_pandas(res, schema=self.schema)

import os

import pyarrow as pa
import pandas as pd
import numpy as np
from pyteomics import mzml
from Bio import Entrez
from Bio import SeqIO
import codecs
import gc

tmt = {
    'TMT10':  ['TMT126', 'TMT127C', 'TMT127N', 'TMT128C', 'TMT128N', 'TMT129C', 'TMT129N', 'TMT130C', 'TMT130N', 'TMT131'],
    'TMT11':  ["TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N", "TMT130C", "TMT131N", "TMT131C"],
    'TMT16':  ["TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N", "TMT130C", "TMT131N", "TMT131C", "TMT132N", "TMT132C", "TMT133N", "TMT133C", "TMT134N"],
    'TMT6':   ["TMT126", "TMT127", "TMT128", "TMT129", "TMT130", "TMT131"]
}


class FeatureConvertor():
    def __init__(self, experiment_type, schema):
        self.experiment_type = experiment_type
        self.schema = schema
        self.mzml_directory = None
        self.email = None
        self.tool = None

    # load csv from mzTab
    def __skip_and_load_csv(self, fle, header, **kwargs):
        fle_len = self.__extract_len(fle, header)
        if os.stat(fle).st_size == 0:
            raise ValueError("File is empty")
        with open(fle) as f:
            pos = 0
            line = f.readline()
            while line.split("\t")[0] != header:
                pos = f.tell()
                line = f.readline()
            f.seek(pos)
            return pd.read_csv(f, nrows=fle_len, **kwargs)

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

    # extract ms runs
    def extract_ms_runs(self, fle):  # DONE implemented.
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

    # extract from mzTab
    def _extract_from_mzTab(self, mzTab_path):
        # load ms_runs
        ms_runs = self.extract_ms_runs(mzTab_path)
        # load PSM
        PSM = self.__skip_and_load_csv(mzTab_path, 'PSH', sep='\t')
        # load PEP
        if self.experiment_type != "lfq":
            PEP = self.__skip_and_load_csv(mzTab_path, 'PEH', sep='\t', usecols=[
                                           "retention_time", 'accession', "opt_global_cv_MS:1000889_peptidoform_sequence", "best_search_engine_score[1]", "charge"])
        else:
            PEP = self.__skip_and_load_csv(mzTab_path, 'PEH', sep='\t', usecols=[
                                           'spectra_ref', 'accession', "opt_global_cv_MS:1000889_peptidoform_sequence", "best_search_engine_score[1]", "charge"])
        peptidoform_sequence = "".join(
            filter(lambda x: x.endswith('peptidoform_sequence'), PSM.columns))
        decoy_peptide = "".join(
            filter(lambda x: x.endswith('decoy_peptide'), PSM.columns))

        # compute count
        spectral_count = pd.DataFrame(PSM.groupby(
            ['sequence', 'charge', peptidoform_sequence]).size()).reset_index()
        spectral_count.columns = ['sequence', 'charge',
                                  peptidoform_sequence, 'spectral_count']
        PSM = PSM.merge(spectral_count, on=[
                        'sequence', 'charge', peptidoform_sequence], how='left')
        # merge PEP
        if self.experiment_type != "lfq":
            PSM['retention_time'] = PSM['retention_time'].round(1)
            PEP['retention_time'] = PEP['retention_time'].round(1)
            PSM = PSM.merge(PEP, on=[
                            'accession', 'charge', 'opt_global_cv_MS:1000889_peptidoform_sequence', 'retention_time'], how='right')
        else:
            PSM['spectra_ref'] = PSM['spectra_ref'].apply(
                lambda x: ms_runs[x.split(":")[0]])
            PEP['spectra_ref'] = PEP['spectra_ref'].apply(
                lambda x: ms_runs[x.split(":")[0]])
            PSM = PSM.merge(PEP, on=[
                            'accession', 'charge', 'opt_global_cv_MS:1000889_peptidoform_sequence', 'spectra_ref'], how='right')
        # extract columns and rename
        PSM = PSM[['sequence',
                   'accession',
                   'unique',
                   'search_engine_score[1]',
                   'modifications',
                   'retention_time',
                   'charge',
                   'exp_mass_to_charge',
                   'calc_mass_to_charge',
                   'start',
                   'end',
                   'opt_global_Posterior_Error_Probability_score',
                   peptidoform_sequence,
                   decoy_peptide,
                   'spectral_count',
                   "best_search_engine_score[1]",
                   'spectra_ref'
                   ]]
        PSM = PSM.rename(columns={
            'accession': 'protein_accessions',
            'start': 'protein_start_positions',
            'end': 'protein_end_positions',
                   peptidoform_sequence: 'peptidoform',
                   'opt_global_Posterior_Error_Probability_score': 'posterior_error_probability',
                   'search_engine_score[1]': 'global_qvalue',
                   decoy_peptide: 'is_decoy',
                   "best_search_engine_score[1]": 'best_id_score',
            'spectra_ref': 'reference_file_name'
        })
        return PSM

    # extract from msstats_in
    def _extract_from_msstats(self, msstats_path):
        msstats_in = pd.read_csv(msstats_path)
        msstats_in['Reference'] = msstats_in['Reference'].apply(
            lambda x: x.split(".")[0])
        if 'FragmentIon' not in msstats_in.columns:
            msstats_in.loc[:, 'FragmentIon'] = 'NA'
        if 'IsotopeLabelType' not in msstats_in.columns:
            msstats_in.loc[:, 'IsotopeLabelType'] = 'L'

        if self.experiment_type != 'lfq':
            msstats_in["Channel"] = msstats_in["Channel"].apply(
                lambda row: tmt[self.experiment_type][row-1])
            msstats_in["RetentionTime"] = msstats_in["RetentionTime"].round(1)
            msstats_in.rename(columns={
                'Charge': 'charge',
                'RetentionTime': 'retention_time'
            }, inplace=True)

        # load tmt columns
        elif "Channel" not in msstats_in.columns:
            msstats_in.loc[:, 'Channel'] = 'label free sample'
            msstats_in.rename(columns={
                'PrecursorCharge': 'charge'
            }, inplace=True)

        msstats_in = msstats_in.rename(columns={
            'ProteinName': 'protein_accessions',
            'PeptideSequence': 'peptidoform',
            'Reference': 'reference_file_name',
            'Channel': 'channel',
            'Run': 'run',
            'BioReplicate': 'biological_replicate',
            'Intensity': 'intensity',
            'FragmentIon': 'fragment_ion',
            'IsotopeLabelType': 'isotope_label_type',
        })

        if self.experiment_type != 'lfq':
            return msstats_in[['protein_accessions', 'peptidoform', 'channel', 'run', 'biological_replicate', 'intensity', 'fragment_ion', 'isotope_label_type', 'charge', 'retention_time']]
        else:
            return msstats_in[['protein_accessions', 'peptidoform', 'reference_file_name', 'channel', 'run', 'biological_replicate', 'intensity', 'fragment_ion', 'isotope_label_type', 'charge']]

    def _extract_from_sdrf(self, sdrf_path):
        sdrf = pd.read_csv(sdrf_path, sep='\t')
        sdrf['comment[data file]'] = sdrf['comment[data file]'].apply(
            lambda x: x.split(".")[0])
        sdrf = sdrf.rename(columns={
            'comment[data file]': 'reference_file_name',
            'source name': 'sample_accession',
            'factor value[organism part]': 'condition',
            'comment[fraction identifier]': 'fraction',
            'comment[label]': 'channel'
        })
        # extract
        if self.experiment_type != 'lfq':
            sdrf = sdrf[['reference_file_name', 'sample_accession',
                         'condition', 'fraction', 'channel']]
            return sdrf
        else:
            sdrf = sdrf[['reference_file_name',
                         'sample_accession', 'condition', 'fraction']]
            return sdrf

    def extract_optional_spectrum_featrue(self, res_df, mzml_directory):

        self.mzml_directory = mzml_directory
        res_df[['mz', 'array_intensity', 'num_peaks']] = res_df[['spectra_ref', 'scan_number']].apply(
            lambda x: self._map_spectrum_mz(x['spectra_ref'], x['scan_number']), axis=1, result_type="expand")

        return res_df

    def _map_spectrum_mz(self, mz_path, scan):
        scan = "controllerType=0 controllerNumber=1 scan=" + str(scan)
        mz_path = self.mzml_directory + '/' + mz_path + '.mzml'
        spectral = mzml.MzML(mz_path).get_by_id(scan)
        mz_array = spectral["m/z array"]
        array_intensity = spectral['intensity array']
        return mz_array, array_intensity, len(mz_array)

    def _query_accession(self, gene_name):
        Entrez.email = self.Email
        Entrez.tool = self.tool
        if gene_name == 'unknow':
            return ['unknow']
        q_name = gene_name + " AND Homo sapiens[porgn:__txid9606]"
        handle = Entrez.esearch(db="nucleotide", term=q_name)
        record = Entrez.read(handle)
        id_list = record["IdList"]
        accessions = self._get_accessions(id_list)
        return accessions

    def _get_accessions(self, id_list):
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
                lambda x: self._query_accession(x))
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

    def __clear_up_memory(self, *dfs):
        for df in dfs:
            del df
            gc.collect()

    def merge_to_table(self, mzTab_path, msstats_path, sdrf_path):

        # load
        psm = self._extract_from_mzTab(mzTab_path)
        msstats = self._extract_from_msstats(msstats_path)
        sdrf = self._extract_from_sdrf(sdrf_path)
        # merge psm and msstats
        if self.experiment_type != 'lfq':
            res = pd.merge(psm, msstats, on=[
                           'protein_accessions', 'peptidoform', "charge", 'retention_time'], how="inner")
            self.__clear_up_memory(psm, msstats)
            res = pd.merge(
                res, sdrf, on=['reference_file_name', 'channel'], how='left')
        else:
            res = pd.merge(psm, msstats, on=[
                           'protein_accessions', 'peptidoform', "charge", 'reference_file_name'], how="inner")
            self.__clear_up_memory(psm, msstats)
            res = pd.merge(res, sdrf, on='reference_file_name', how='left')

        # optional fields
        res.loc[:, 'gene_accessions'] = "not"
        res.loc[:, 'gene_names'] = "not"
        res.loc[:, 'id_scores'] = "not"
        res.loc[:, 'mz'] = 0.0
        res.loc[:, 'intensity_array'] = 0.0
        res.loc[:, 'num_peaks'] = 0

        res['gene_accessions'] = res['gene_accessions'].apply(lambda x: [x])
        res['gene_names'] = res['gene_names'].apply(lambda x: [x])
        res['id_scores'] = res['id_scores'].apply(lambda x: [x])
        res['mz'] = res['mz'].apply(lambda x: [x])
        res['intensity_array'] = res['intensity_array'].apply(lambda x: [x])

        return res

    def convert_to_parquet(self, res):
        res['sequence'] = res['sequence'].fillna("").astype(str)
        res['protein_accessions'] = res['protein_accessions'].str.split('|')
        res['unique'] = res['unique'].fillna(0).astype(int)
        res['modifications'] = res['modifications'].fillna("").str.split(',')
        res['retention_time'] = res['retention_time'].fillna(0.0)
        res['charge'] = res['charge'].fillna(0).astype(int)
        res['exp_mass_to_charge'] = res['exp_mass_to_charge']
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
        res['scan_number'].astype(str)
        res['consensus_support'] = res['consensus_support'].fillna(0.0)

        return pa.Table.from_pandas(res, schema=self.schema)

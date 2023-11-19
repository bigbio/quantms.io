import os
from unittest import TestCase

import click
import pandas as pd
from click.testing import CliRunner

from quantms_io.commands.diann_convert_command import diann_convert_to_parquet
from quantms_io.core.diann_convert import get_exp_design_dfs, find_modification, generate_scan_number, \
    handle_protein_map, DiaNNConvert
import time

from quantms_io.core.project import create_uuid_filename


class Test(TestCase):


    def test_all_required_arguments(self):
        st = time.time()

        # Initialize input arguments
        report_path = "/Users/yperez/work/quantms-data/PXD037340.2/diann_report.tsv"
        design_file = "/Users/yperez/work/quantms-data/PXD037340.2/PXD037340-DIA.sdrf_openms_design.tsv"
        fasta_path = "/Users/yperez/work/multiomics-configs/databases/Homo-sapiens-uniprot-reviewed-contaminants-202210.fasta"
        modifications = ["Carbamidomethyl (C)", "Oxidation (M)"]
        pg_path = "/Users/yperez/work/quantms-data/PXD037340.2/diann_report.pg_matrix.tsv"
        pr_path = "/Users/yperez/work/quantms-data/PXD037340.2/diann_report.pr_matrix.tsv"
        qvalue_threshold = 0.01
        mzml_info_folder = "/Users/yperez/work/quantms-data/PXD037340.2/mzmlstatistics"
        sdrf_path = "/Users/yperez/work/quantms-data/PXD037340.2/PXD037340-DIA.sdrf.tsv"
        output_folder = "/Users/yperez/work/quantms-data/PXD037340.2/"
        output_prefix_file = "PXD037340.2"
        chunksize = 100000

        # Invoke the function
        # Assert that the feature and psm files are generated in the specified output folder
        feature_output_path = output_folder + "/" + create_uuid_filename(output_prefix_file, '.feature.parquet')
        psm_output_path = output_folder + "/" + create_uuid_filename(output_prefix_file, '.psm.parquet')

        DiaNN = DiaNNConvert()

        DiaNN.generate_feature_and_psm_file(report_path=report_path, design_file=design_file, fasta_path=fasta_path,
                                            modifications=modifications, pg_path=pg_path, pr_path=pr_path,
                                            qvalue_threshold=qvalue_threshold, mzml_info_folder=mzml_info_folder,
                                            sdrf_path=sdrf_path, output_path=feature_output_path,
                                            psm_output_path=psm_output_path, chunksize=chunksize)
        et = time.time()

        # get the execution time
        elapsed_time = et - st
        print('Execution time:', elapsed_time, 'seconds')

    def test_read_experimental_design_file(self):
        # Measure the time to read the file
        # get the start time
        st = time.time()

        exp_design_file = "data/PXD030304.sdrf_openms_design.tsv"

        # Act
        result = get_exp_design_dfs(exp_design_file)

        # Assert
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], pd.DataFrame)
        assert isinstance(result[1], pd.DataFrame)

        # get the end time
        et = time.time()

        # get the execution time
        elapsed_time = et - st
        print('Execution time:', elapsed_time, 'seconds')

        # Check that the file name of the firt entry in fraction table is 190204_2300_0054N_008SA_M06_S_1
        assert result[1].iloc[0]['run'] == '190204_2300_0054N_008SA_M06_S_1'

    def test_single_modification(self):
        st = time.time()
        assert find_modification("PEPM(UNIMOD:35)IDE") == "4-UNIMOD:35"
        et = time.time()

        # get the execution time
        elapsed_time = et - st
        print('Execution time:', elapsed_time, 'seconds')

    def test_returns_scan_number_when_scan_in_spectra_ref(self):
        spectra_ref = 'scan=123'
        assert generate_scan_number(spectra_ref) == '123'

    def test_returns_value_if_key_exists(self):
        protein_map = {"key1": "value1", "key2": "value2"}
        key = "key1"
        result = handle_protein_map(protein_map, key)
        assert result == "value1"

    def test_update_dict_with_new_key_value_pairs(self):
        map_dict = {'A': 1.2, 'B': 3.2}
        temporary_dict = {'A': 2.8, 'C': 3.3}

        converter = DiaNNConvert()
        updated_dict = converter._DiaNNConvert__update_dict(map_dict, temporary_dict, 0)

        assert updated_dict == {'A': 1.2, 'B': 3.2, 'C': 3.3}


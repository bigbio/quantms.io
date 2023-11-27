from unittest import TestCase
import pandas as pd

from quantms_io.core.diann_convert import get_exp_design_dfs, find_modification, DiaNNConvert
import time

from quantms_io.core.project import create_uuid_filename


class Test(TestCase):


    def test_all_required_arguments(self):
        st = time.time()

        # Initialize input arguments
        report_path = "/Users/yperez/work/quantms-data/PXD037340.2/diann_report.tsv"
        design_file = "/Users/yperez/work/quantms-data/PXD037340.2/PXD037340-DIA.sdrf_openms_design.tsv"
        sdrf_file   = "/Users/yperez/work/quantms-data/PXD037340.2/PXD037340-DIA.sdrf.tsv"
        modifications = ["Carbamidomethyl (C)", "Oxidation (M)"]
        qvalue_threshold = 0.01
        mzml_info_folder = "/Users/yperez/work/quantms-data/PXD037340.2/mzmlstatistics"
        output_folder = "/Users/yperez/work/quantms-data/PXD037340.2/"
        output_prefix_file = "PXD037340.2"


        # Invoke the function
        # Assert that the feature and psm files are generated in the specified output folder
        feature_output_path = output_folder + "/" + create_uuid_filename(output_prefix_file, '.feature.parquet')
        psm_output_path = output_folder + "/" + create_uuid_filename(output_prefix_file, '.psm.parquet')

        DiaNN = DiaNNConvert()

        DiaNN.generate_psm_and_feature_file(report_path=report_path, design_file=design_file,
                                            modifications=modifications, qvalue_threshold=qvalue_threshold,
                                            feature_output_path=feature_output_path,
                                            psm_output_path=psm_output_path, sdrf_path=sdrf_file,
                                            mzml_info_folder=mzml_info_folder, file_num=10,
                                            duckdb_max_memory="4GB", duckdb_threads=4)
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


from unittest import TestCase

import pandas as pd
from quantms_io.core.diann_convert import get_exp_design_dfs, find_modification, generate_scan_number, \
    handle_protein_map
import time

class Test(TestCase):

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
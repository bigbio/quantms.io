import datetime
import os
import time
from unittest import TestCase

from quantms_io.core.mztab import MztabHandler
from quantms_io.core.psm import PSMHandler


class TestMztabHandler(TestCase):
    def test_load_mztab_file(self):
        # get the current time of the machine
        time = datetime.datetime.now()
        mztab_handler = MztabHandler(
            mztab_file="data/raw_ae_example/PXD009219.sdrf_openms_design_openms.mzTab",
            use_cache=True,
        )
        mztab_handler.load_mztab_file(use_cache=True)
        mztab_handler.print_mztab_stats()
        mztab_handler.close()
        print("Time to load the mztab file: {}".format(datetime.datetime.now() - time))

    def test_valid_mztab_conversion(self):
        # Initialize

        st = time.time()
        mztab_path = "/Users/yperez/work/quantms.io/python/quantmsio/quantms_io/data/raw_ae_example/PXD009219.sdrf_openms_design_openms.mzTab"
        parquet_path = "/Users/yperez/work/quantms.io/python/quantmsio/quantms_io/data/PXD009219.psm.parquet"
        verbose = False
        batch_size = 100000

        psm_handler = PSMHandler()

        # Invoke
        psm_handler.convert_mztab_to_psm(mztab_path=mztab_path, parquet_path = parquet_path, verbose=verbose, batch_size=batch_size)

        # Assert
        assert os.path.exists(parquet_path)
        assert os.path.isfile(parquet_path)

        et = time.time()

        # get the execution time
        elapsed_time = et - st
        print('Execution time:', elapsed_time, 'seconds')
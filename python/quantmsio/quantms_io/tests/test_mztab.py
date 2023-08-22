import datetime
from unittest import TestCase

from quantms_io.core.mztab import MztabHandler


class TestMztabHandler(TestCase):
    def test_load_mztab_file(self):
        # get the current time of the machine
        time = datetime.datetime.now()
        # mztab_handler = MztabHandler(mztab_file="/Users/yperez/Downloads/PXD002137.sdrf_openms_design_openms.mzTab", use_cache=True)
        mztab_handler = MztabHandler(mztab_file="data/raw_ae_example/PXD009219.sdrf_openms_design_openms.mzTab", use_cache=True)
        mztab_handler.load_mztab_file(use_cache=True)
        mztab_handler.print_mztab_stats()
        mztab_handler.close()
        print("Time to load the mztab file: {}".format(datetime.datetime.now() - time))
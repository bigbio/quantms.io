from unittest import TestCase

from quantms_io.core.psm import PSMHandler


class TestPSMHandler(TestCase):
    def test_convert_mztab_to_feature(self):
        psm_handler = PSMHandler()
        psm_handler.convert_mztab_to_feature(mztab_path="data/raw_ae_example/PXD009219.sdrf_openms_design_openms.mzTab", parquet_path="output.parquet")

from .common import datafile
from unittest import TestCase
from quantmsio.core.psm import Psm


class TestPSMHandler(TestCase):

    def test_convert_mztab_to_feature(self):
        mztab_path = datafile("DDA-lfq/PXD040438.mzTab")
        psm = Psm(mztab_path)
        for _ in psm.generate_report():
            print("ok")

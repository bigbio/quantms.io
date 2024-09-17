from .common import datafile
from unittest import TestCase
from quantmsio.temp_core.psm import PsmInMemory
class TestPSMHandler(TestCase):

    def test_convert_mztab_to_feature(self):
        mztab_path = datafile("DDA-plex/MSV000079033.mzTab")
        psm = PsmInMemory(mztab_path)
        for _ in psm.generate_report():
            print("ok")
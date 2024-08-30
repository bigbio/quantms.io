from unittest import TestCase

from quantmsio.core.psm import PSMHandler
from quantmsio.core.psm_in_memory import PsmInMemory
from .common import datafile


class TestPSMHandler(TestCase):

    def test_convert_mztab_to_feature(self):
        mztab_path = datafile("DDA-plex/MSV000079033.mzTab")

        p = PSMHandler()
        psm = PsmInMemory(p.schema)
        for _ in psm.generate_psm_parquet(mztab_path):
            print("ok")

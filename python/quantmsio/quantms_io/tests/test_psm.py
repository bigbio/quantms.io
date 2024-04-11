from unittest import TestCase

from core.psm import PSMHandler
from core.psm_in_memory import PsmInMemory


class TestPSMHandler(TestCase):

    def test_convert_mztab_to_feature(self):
        mztab_path = __package__ + "/examples/DDA-plex/MSV000079033.mzTab"

        p = PSMHandler()
        psm = PsmInMemory(p.schema)
        for _ in psm.generate_psm_parquet(mztab_path):
            print("ok")

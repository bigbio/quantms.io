from unittest import TestCase

from ddt import data
from ddt import ddt

from quantms_io.core.project import create_uuid_filename
from quantms_io.core.psm import PSMHandler
from quantms_io.core.psm_in_memory import PsmInMemory

#@ddt
class TestPSMHandler(TestCase):
    #global test_data
    
    #@data(*test_data)
    def test_convert_mztab_to_feature(self):
        mztab_path = __package__ + "/examples/DDA-plex/MSV000079033.mzTab"

        p = PSMHandler()
        Psm = PsmInMemory(p.schema)
        for _ in Psm.generate_psm_parquet(mztab_path):
            print('ok')

from unittest import TestCase

from quantms_io.core.psm import PSMHandler
from quantms_io.core.project import create_uuid_filename
from ddt import data,ddt
@ddt
class TestPSMHandler(TestCase):
    global test_datas
    test_datas = [
        ("/examples/DDA-lfq/PXD040438.mzTab","/examples/output/DDA-lfq/","PXD040438"), 
        ("/examples/DDA-plex/MSV000079033.mzTab","/examples/output/DDA-plex/","MSV000079033"),
    ]

    @data(*test_datas)
    def test_convert_mztab_to_feature(self,test_data):
        mztab_path = __package__ + test_data[0]
        output_path = __package__ + test_data[1] + create_uuid_filename(test_data[2],'.psm.parquet')

        psm_handler = PSMHandler()
        psm_handler.convert_mztab_to_psm(
            mztab_path=mztab_path,
            parquet_path=output_path,
        )

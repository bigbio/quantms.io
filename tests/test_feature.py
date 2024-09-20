from .common import datafile
from unittest import TestCase
from quantmsio.core.feature import Feature
from ddt import data
from ddt import ddt
@ddt
class TestFeatureHandler(TestCase):
    global test_datas
    test_datas = [
        (
            "DDA-plex/MSV000079033.mzTab",
            "DDA-plex/MSV000079033_msstats_in.csv",
            "DDA-plex/MSV000079033-Blood-Plasma-iTRAQ.sdrf.tsv",
            "ITRAQ4",
        ),
    ]

    @data(*test_datas)
    def test_transform_msstats(self, test_data):
        mztab_file = datafile(test_data[0])
        msstats_file = datafile(test_data[1])
        sdrf_file = datafile(test_data[2])
        F = Feature(mztab_file,sdrf_file,msstats_file)
        F.transform_msstats_in()
    
    @data(*test_datas)
    def test_extract_psm_msg(self, test_data):
        mztab_file = datafile(test_data[0])
        msstats_file = datafile(test_data[1])
        sdrf_file = datafile(test_data[2])
        F = Feature(mztab_file,sdrf_file,msstats_file)
        F.extract_psm_msg()
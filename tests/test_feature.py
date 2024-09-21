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
        for msstats in F.transform_msstats_in():
            print('ok')
    
    @data(*test_datas)
    def test_extract_psm_msg(self, test_data):
        mztab_file = datafile(test_data[0])
        msstats_file = datafile(test_data[1])
        sdrf_file = datafile(test_data[2])
        F = Feature(mztab_file,sdrf_file,msstats_file)
        F.extract_psm_msg()
    
    @data(*test_datas)
    def test_extract_sdrf(self, test_data):
        mztab_file = datafile(test_data[0])
        msstats_file = datafile(test_data[1])
        sdrf_file = datafile(test_data[2])
        F = Feature(mztab_file,sdrf_file,msstats_file)
        F.transform_sdrf()
    
    @data(*test_datas)
    def test_merge_msstats_and_sdrf(self, test_data):
        mztab_file = datafile(test_data[0])
        msstats_file = datafile(test_data[1])
        sdrf_file = datafile(test_data[2])
        F = Feature(mztab_file,sdrf_file,msstats_file)
        for msstats in F.transform_msstats_in():
            F.merge_msstats_and_sdrf(msstats)
    
    @data(*test_datas)
    def test_merge_msstats_and_psm(self, test_data):
        mztab_file = datafile(test_data[0])
        msstats_file = datafile(test_data[1])
        sdrf_file = datafile(test_data[2])
        F = Feature(mztab_file,sdrf_file,msstats_file)
        map_dict = F.extract_psm_msg()
        for msstats in F.transform_msstats_in():
            F.merge_msstats_and_psm(msstats,map_dict)
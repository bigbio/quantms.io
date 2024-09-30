from .common import datafile
from unittest import TestCase
from quantmsio.core.diann import DiaNNConvert
from ddt import data
from ddt import ddt
@ddt
class TestFeatureHandler(TestCase):
    global test_datas
    test_datas = [
        (
            "DIANN/diann_report.tsv",
            "DIANN/PXD019909-DIA.sdrf_openms_design.tsv",
            "DIANN/PXD019909-DIA.sdrf.tsv",
            "DIANN/mzml"
        ),
    ]

    @data(*test_datas)
    def test_transform_msstats(self, test_data):
        report_file = datafile(test_data[0])
        design_file = datafile(test_data[1])
        sdrf_file = datafile(test_data[2])
        mzml = datafile(test_data[3])
        D = DiaNNConvert(report_file,sdrf_file)
        for msstats in D.main_report_df(0.05,mzml,2):
            print('ok')
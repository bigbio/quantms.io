from .common import datafile
from unittest import TestCase
from quantmsio.core.diann import DiaNNConvert
from quantmsio.core.feature import Feature
from ddt import data
from ddt import ddt


@ddt
class TestFeatureHandler(TestCase):
    test_datas = [
        (
            "DIANN/diann_report.tsv",
            "DIANN/PXD019909-DIA.sdrf.tsv",
            "DIANN/mzml",
        ),
    ]

    @data(*test_datas)
    def test_transform_feature(self, test_data):
        report_file = datafile(test_data[0])
        sdrf_file = datafile(test_data[1])
        mzml = datafile(test_data[2])
        D = DiaNNConvert(report_file, sdrf_file)
        for report in D.main_report_df(0.05, mzml, 2):
            D.add_additional_msg(report)
            Feature.convert_to_parquet_format(report)
            Feature.transform_feature(report)

    @data(*test_datas)
    def test_transform_features(self, test_data):
        report_file = datafile(test_data[0])
        sdrf_file = datafile(test_data[1])
        mzml = datafile(test_data[2])
        D = DiaNNConvert(report_file, sdrf_file)
        for report in D.main_report_df(0.05, mzml, 2):
            D.add_additional_msg(report)
            Feature.convert_to_parquet_format(report)
            for _, df in Feature.slice(report, ["reference_file_name", "precursor_charge"]):
                Feature.transform_feature(df)

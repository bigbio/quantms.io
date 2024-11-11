from quantmsio.core.feature import Feature
from quantmsio.core.maxquant import MaxQuant
from quantmsio.core.sdrf import SDRFHandler
from .common import datafile
from unittest import TestCase
from quantmsio.core.psm import Psm
from ddt import data
from ddt import ddt


@ddt
class TestFeatureHandler(TestCase):
    test_datas = [
        (
            "maxquant/msms.txt",
            "maxquant/evidence.txt",
            "maxquant/sdrf.tsv",
        ),
    ]

    @data(*test_datas)
    def test_transform_psm(self, test_data):
        msms_file = datafile(test_data[0])
        M = MaxQuant()
        for df in M.iter_batch(msms_file, "psm", chunksize=500000):
            M.transform_psm(df)
            Psm.convert_to_parquet_format(df)
            Psm.transform_parquet(df)

    @data(*test_datas)
    def test_transform_feature(self, test_data):
        evidence_file = datafile(test_data[1])
        sdrf_file = datafile(test_data[2])
        Sdrf = SDRFHandler(sdrf_file)
        M = MaxQuant()
        M.experiment_type = Sdrf.get_experiment_type_from_sdrf()
        M._sample_map = Sdrf.get_sample_map_run()
        for df in M.iter_batch(evidence_file, chunksize=500000):
            M.transform_feature(df)
            Feature.convert_to_parquet_format(df)
            Feature.transform_feature(df)

    @data(*test_datas)
    def test_transform_features(self, test_data):
        evidence_file = datafile(test_data[1])
        sdrf_file = datafile(test_data[2])
        Sdrf = SDRFHandler(sdrf_file)
        M = MaxQuant()
        M.experiment_type = Sdrf.get_experiment_type_from_sdrf()
        M._sample_map = Sdrf.get_sample_map_run()
        for report in M.iter_batch(evidence_file, chunksize=500000):
            M.transform_feature(report)
            Feature.convert_to_parquet_format(report)
            for _, df in Feature.slice(report, ["reference_file_name","precursor_charge"]):
                Feature.transform_feature(df)
            
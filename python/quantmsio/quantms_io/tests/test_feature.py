import random
import unittest
from unittest import TestCase

from ddt import data
from ddt import ddt

from quantms_io.core.feature import FeatureHandler
from quantms_io.core.project import create_uuid_filename


@ddt
class TestFeatureHandler(TestCase):
    global test_datas
    test_datas = [
        (
            "/examples/DDA-lfq/PXD040438.mzTab",
            "/examples/DDA-lfq/PXD040438_msstats_in.csv",
            "/examples/DDA-lfq/PXD040438.sdrf.tsv",
            "/examples/output/DDA-lfq/",
            "PXD040438",
        ),
        (
            "/examples/DDA-plex/MSV000079033.mzTab",
            "/examples/DDA-plex/MSV000079033_msstats_in.csv",
            "/examples/DDA-plex/MSV000079033-Blood-Plasma-iTRAQ.sdrf.tsv",
            "/examples/output/DDA-plex/",
            "MSV000079033",
        ),
    ]

    @data(*test_datas)
    def test_convert_mztab_msstats_to_feature(self, test_data):
        mztab_file = __package__ + test_data[0]
        msstats_file = __package__ + test_data[1]
        sdrf_file = __package__ + test_data[2]

        feature_manager = FeatureHandler()
        feature_manager.parquet_path = (
            __package__ + test_data[3] + create_uuid_filename(test_data[4], ".feature.parquet")
        )
        feature_manager.convert_mztab_msstats_to_feature(
            mztab_file=mztab_file,
            msstats_file=msstats_file,
            sdrf_file=sdrf_file,
            use_cache=False,
        )

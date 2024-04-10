from unittest import TestCase

from quantms_io.core.json import JsonConverter


class TestJson(TestCase):

    def test_project_file_to_json(self):
        feature_path = ["/examples/JSON/PXD040438.feature.parquet", "/examples/output/JSON/PXD040438.featrue.json"]
        ae_path = ["/examples/JSON/PXD016999.absolute.tsv", "/examples/output/JSON/PXD016999.ae.json"]
        sdrf_path = [
            "/examples/JSON/MSV000079033-Blood-Plasma-iTRAQ.sdrf.tsv",
            "/examples/output/JSON/MSV000079033.sdrf.json",
        ]
        feature_path = [__package__ + path for path in feature_path]
        ae_path = [__package__ + path for path in ae_path]
        sdrf_path = [__package__ + path for path in sdrf_path]
        converter = JsonConverter()

        converter.psm_to_json(feature_path[0], feature_path[1])
        converter.convert_tsv_to_json(ae_path[0], ae_path[1])
        converter.sdrf_to_json(sdrf_path[0], sdrf_path[1])
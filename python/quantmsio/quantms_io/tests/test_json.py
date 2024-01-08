from unittest import TestCase
from quantms_io.core.json import JsonConverter

class TestJson(TestCase):

    def test_project_file_to_json(self):
        feature_path = ['data\PXD014414\PXD014414-943a8f02-0527-4528-b1a3-b96de99ebe75.featrue.parquet','data\json\PXD014414.featrue.json']
        psm_path = ['data\PXD014414\PXD014414-f4fb88f6-0a45-451d-a8a6-b6d58fb83670.psm.parquet','data\json\PXD014414.psm.json']
        de_path = ['data\PXD014414\PXD014414-48fc92e9-ffd2-4335-a6f5-ce7d5278b5ff.differential.tsv','data\json\PXD014414.de.json']
        ae_path = ['data\raw_ae_example\quantmsio\PXD016999.1-c6fa5b54-1257-4b40-a42c-8754d01a657b.absolute.tsv','data\json\PXD016999.1.ae.json']
        converter = JsonConverter()

        converter.psm_to_json(feature_path[0], feature_path[1])
        converter.feature_to_json(psm_path[0], psm_path[1])
        converter.convert_tsv_to_json(de_path[0],de_path[1])
        converter.convert_tsv_to_json(ae_path[0],ae_path[1])

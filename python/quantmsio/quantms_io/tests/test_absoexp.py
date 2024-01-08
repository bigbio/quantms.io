from unittest import TestCase

from quantms_io.core.ae import AbsoluteExpressionHander


class TestAEHandler(TestCase):

    def test_write_an_example_ae(self):
        ae_handler = AbsoluteExpressionHander()
        ae_handler.load_project_file("data/raw_ae_example/quantmsio/project.json")
        ae_handler.load_ibaq_file("data\raw_ae_example\PXD016999.1-ibaq.tsv")
        ae_handler.load_sdrf_file("data\raw_ae_example\PXD016999-first-instrument.sdrf.tsv")
        ae_handler.convert_ibaq_to_quantms(
            output_folder="data/raw_ae_example/quantmsio/",
            output_file_prefix="PXD016999.1",
            delete_existing=True,
        )
        

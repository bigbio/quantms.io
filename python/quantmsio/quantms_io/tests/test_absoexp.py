from unittest import TestCase

from core.ae import AbsoluteExpressionHander


class TestAEHandler(TestCase):

    def test_ae(self):
        test_data = (
            "/examples/AE/project.json",
            "/examples/AE/PXD016999.1-ibaq.tsv",
            "/examples/AE/PXD016999-first-instrument.sdrf.tsv",
        )

        project_path = __package__ + test_data[0]
        ibaq_path = __package__ + test_data[1]
        sdrf_path = __package__ + test_data[2]
        ae_handler = AbsoluteExpressionHander()
        ae_handler.load_project_file(project_path)
        ae_handler.load_ibaq_file(ibaq_path)
        ae_handler.load_sdrf_file(sdrf_path)

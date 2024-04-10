from unittest import TestCase

from quantms_io.core.ae import AbsoluteExpressionHander


class TestAEHandler(TestCase):

    def test_write_an_example_ae(self):
        test_data = (
            "/examples/AE/project.json",
            "/examples/AE/PXD016999.1-ibaq.tsv",
            "/examples/AE/PXD016999-first-instrument.sdrf.tsv",
        )

        project_path = __package__ + test_data[0]
        ibaq_path = __package__ + test_data[1]
        sdrf_path = __package__ + test_data[2]

        output_folder = __package__ + "/examples/output/AE/"

        output_file_prefix = "PXD016999.1"
        ae_handler = AbsoluteExpressionHander()
        ae_handler.load_project_file(project_path)
        ae_handler.load_ibaq_file(ibaq_path)
        ae_handler.load_sdrf_file(sdrf_path)
        ae_handler.convert_ibaq_to_quantms(
            output_folder=output_folder,
            output_file_prefix=output_file_prefix,
            delete_existing=False,
        )

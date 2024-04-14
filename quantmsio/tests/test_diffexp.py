from unittest import TestCase

from core.de import DifferentialExpressionHandler


class TestDEHandler(TestCase):

    def test_write_an_example_de(self):
        test_data = (
            "/examples/DE/project.json",
            "/examples/DE/PXD033169.sdrf_openms_design_msstats_in_comparisons.csv",
            "/examples/DE/PXD033169.sdrf.tsv",
        )

        project_path = __package__ + test_data[0]
        de_path = __package__ + test_data[1]
        sdrf_path = __package__ + test_data[2]
        de_handler = DifferentialExpressionHandler()
        de_handler.load_project_file(project_path)
        de_handler.load_msstats_file(de_path)
        de_handler.load_sdrf_file(sdrf_path)
        de_handler.set_fdr_threshold(fdr_threshold=0.05)

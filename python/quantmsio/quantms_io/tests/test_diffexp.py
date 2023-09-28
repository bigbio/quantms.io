from unittest import TestCase

from quantms_io.core.de import DifferentialExpressionHandler


class TestDEHandler(TestCase):

    def test_write_an_example_de(self):
        de_handler = DifferentialExpressionHandler()
        de_handler.load_project_file("data/raw_de_example/quantmsio/project.json")
        de_handler.load_msstats_file("data/raw_de_example/PXD033169.sdrf_openms_design_msstats_in_comparisons.csv")
        de_handler.load_sdrf_file("data/raw_de_example/PXD033169.sdrf.tsv")
        de_handler.set_fdr_threshold(fdr_threshold=0.05)
        de_handler.convert_msstats_to_quantms(
            output_folder="data/raw_de_example/quantmsio/",
            output_file_prefix="PXD033169",
            delete_existing=True,
        )
        de_handler.update_project_file("data/raw_de_example/quantmsio/project.json")




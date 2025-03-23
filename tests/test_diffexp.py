from pathlib import Path

from quantmsio.core.de import DifferentialExpressionHandler

TEST_DATA_ROOT = Path(__file__).parent / "examples"


def test_differential_expression():
    """Test differential expression handler."""
    # Resolve file paths
    test_data = (
        TEST_DATA_ROOT / "DE/project.json",
        TEST_DATA_ROOT / "DE/PXD033169.sdrf_openms_design_msstats_in_comparisons.csv",
        TEST_DATA_ROOT / "DE/PXD033169.sdrf.tsv",
    )

    project_path = test_data[0]
    de_path = test_data[1]
    sdrf_path = test_data[2]
    de_handler = DifferentialExpressionHandler()
    de_handler.load_project_file(project_path)
    de_handler.load_msstats_file(de_path)
    de_handler.load_sdrf_file(sdrf_path)
    de_handler.set_fdr_threshold(fdr_threshold=0.05)

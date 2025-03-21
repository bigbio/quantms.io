import pytest
from quantmsio.core.de import DifferentialExpressionHandler
from .common import datafile


def test_differential_expression():
    """Test differential expression handler."""
    # Resolve file paths
    test_data = (
        "DE/project.json",
        "DE/PXD033169.sdrf_openms_design_msstats_in_comparisons.csv",
        "DE/PXD033169.sdrf.tsv",
    )

    project_path = datafile(test_data[0])
    de_path = datafile(test_data[1])
    sdrf_path = datafile(test_data[2])

    # Initialize DifferentialExpressionHandler
    de_handler = DifferentialExpressionHandler()

    # Load project file
    de_handler.load_project_file(project_path)
    assert de_handler.project is not None
    assert "project_accession" in de_handler.project.project_info

    # Load msstats file
    de_handler.load_msstats_file(de_path)
    assert de_handler.de_df is not None
    assert len(de_handler.de_df) > 0
    assert "Protein" in de_handler.de_df.columns
    assert "log2FC" in de_handler.de_df.columns
    assert "pvalue" in de_handler.de_df.columns

    # Load sdrf file
    de_handler.load_sdrf_file(sdrf_path)
    assert de_handler.sdrf is not None
    assert hasattr(de_handler.sdrf, "get_sample_map")

    # Set FDR threshold
    de_handler.set_fdr_threshold(fdr_threshold=0.05)
    assert de_handler.fdr_threshold == 0.05

    # Test generating differential expression
    de_df = de_handler.generate_differential_expression()
    assert de_df is not None
    assert len(de_df) > 0
    assert "Protein" in de_df.columns
    assert "log2FC" in de_df.columns
    assert "pvalue" in de_df.columns
    assert "significant" in de_df.columns

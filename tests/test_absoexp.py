from quantmsio.core.ae import AbsoluteExpressionHander
from .common import datafile


def test_ae():
    """Test absolute expression handler."""
    # Resolve file paths
    project_path, ibaq_path, sdrf_path = map(
        datafile,
        (
            "AE/project.json",
            "AE/PXD016999.1-ibaq.tsv",
            "AE/PXD016999-first-instrument.sdrf.tsv",
        ),
    )

    # Initialize AbsoluteExpressionHander
    ae_handler = AbsoluteExpressionHander()

    # Load a project file
    ae_handler.load_project_file(project_path)
    assert ae_handler.project_manager is not None
    assert "project_accession" in ae_handler.project_manager.project.project_info

    # Load ibaq file
    ae_handler.load_ibaq_file(ibaq_path)
    assert ae_handler.ibaq_df is not None
    assert len(ae_handler.ibaq_df) > 0
    assert "protein" in ae_handler.ibaq_df.columns
    assert "ibaq" in ae_handler.ibaq_df.columns

    # Load sdrf file
    ae_handler.load_sdrf_file(sdrf_path)
    assert ae_handler.sdrf_manager is not None
    assert hasattr(ae_handler.sdrf_manager, "get_sample_map")

    # Test generating absolute expression
    assert ae_handler.ibaq_df is not None
    assert len(ae_handler.ibaq_df) > 0
    assert "protein" in ae_handler.ibaq_df

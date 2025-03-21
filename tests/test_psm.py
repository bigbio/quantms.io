import pytest
from .common import datafile
from quantmsio.core.psm import Psm


def test_convert_mztab_to_feature():
    """Test converting mzTab to feature."""
    # Resolve file path
    mztab_path = datafile("DDA-lfq/PXD040438.mzTab")

    # Initialize Psm
    psm = Psm(mztab_path)

    # Generate report
    count = 0
    for report in psm.generate_report():
        # Add assertions to verify the result
        assert report is not None
        assert len(report) > 0
        assert "peptidoform" in report.columns
        count += 1

    # Ensure we got at least one result
    assert count > 0


def test_extract_from_pep():
    """Test extracting data from PEP section."""
    # Resolve file path
    mztab_path = datafile("DDA-lfq/PXD040438.mzTab")

    # Initialize Psm
    psm = Psm(mztab_path)

    # Extract from PEP
    pep_dict = psm.extract_from_pep()

    # Add assertions to verify the result
    assert pep_dict is not None
    assert isinstance(pep_dict, dict)
    assert len(pep_dict) > 0


def test_iter_psm_table():
    """Test iterating through PSM table."""
    # Resolve file path
    mztab_path = datafile("DDA-lfq/PXD040438.mzTab")

    # Initialize Psm
    psm = Psm(mztab_path)

    # Iterate through PSM table
    count = 0
    for psm_df in psm.iter_psm_table():
        # Add assertions to verify the result
        assert psm_df is not None
        assert len(psm_df) > 0
        assert "peptidoform" in psm_df.columns
        count += 1

    # Ensure we got at least one result
    assert count > 0


def test_get_mods_map():
    """Test getting modifications map."""
    # Resolve file path
    mztab_path = datafile("DDA-lfq/PXD040438.mzTab")

    # Initialize Psm
    psm = Psm(mztab_path)

    # Get modifications map
    mods_map = psm.get_mods_map()

    # Add assertions to verify the result
    assert mods_map is not None
    assert isinstance(mods_map, dict)
    assert len(mods_map) > 0

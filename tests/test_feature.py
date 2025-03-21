import pytest
from .common import datafile
from quantmsio.core.feature import Feature

test_data = [
    (
        "DDA-lfq/PXD040438.mzTab",
        "DDA-lfq/PXD040438_msstats_in.csv",
        "DDA-lfq/PXD040438.sdrf.tsv",
    ),
]


@pytest.mark.parametrize("mztab_path,msstats_path,sdrf_path", test_data)
def test_transform_msstats(mztab_path, msstats_path, sdrf_path):
    """Test transforming msstats data."""
    # Resolve file paths
    mztab_file = datafile(mztab_path)
    msstats_file = datafile(msstats_path)
    sdrf_file = datafile(sdrf_path)

    # Initialize Feature
    feature = Feature(mztab_file, sdrf_file, msstats_file)

    # Process msstats data
    count = 0
    for msstats in feature.transform_msstats_in():
        # Add assertions to verify the result
        assert msstats is not None
        assert len(msstats) > 0
        assert "peptidoform" in msstats.columns
        count += 1

    # Ensure we got at least one result
    assert count > 0


@pytest.mark.parametrize("mztab_path,msstats_path,sdrf_path", test_data)
def test_extract_psm_msg(mztab_path, msstats_path, sdrf_path):
    """Test extracting PSM messages."""
    # Resolve file paths
    mztab_file = datafile(mztab_path)
    msstats_file = datafile(msstats_path)
    sdrf_file = datafile(sdrf_path)

    # Initialize Feature
    feature = Feature(mztab_file, sdrf_file, msstats_file)

    # Extract PSM messages
    map_dict, pep_dict = feature.extract_psm_msg()

    # Add assertions to verify the result
    assert map_dict is not None
    assert isinstance(map_dict, dict)
    assert len(map_dict) > 0

    assert pep_dict is not None
    assert isinstance(pep_dict, dict)
    assert len(pep_dict) > 0


@pytest.mark.parametrize("mztab_path,msstats_path,sdrf_path", test_data)
def test_generate_feature(mztab_path, msstats_path, sdrf_path):
    """Test generating features."""
    # Resolve file paths
    mztab_file = datafile(mztab_path)
    msstats_file = datafile(msstats_path)
    sdrf_file = datafile(sdrf_path)

    # Initialize Feature
    feature = Feature(mztab_file, sdrf_file, msstats_file)

    # Generate features
    count = 0
    for feature_table in feature.generate_feature():
        # Add assertions to verify the result
        assert feature_table is not None
        assert len(feature_table) > 0
        assert "peptidoform" in feature_table.column_names
        count += 1

    # Ensure we got at least one result
    assert count > 0

from pathlib import Path
from quantmsio.core.feature import Feature

TEST_DATA_ROOT = Path(__file__).parent / "examples"

test_data = (
    TEST_DATA_ROOT / "DDA-lfq/PXD040438.mzTab",
    TEST_DATA_ROOT / "DDA-lfq/PXD040438_msstats_in.csv",
    TEST_DATA_ROOT / "DDA-lfq/PXD040438.sdrf.tsv",
)


def test_transform_msstats():
    """Test transforming msstats data."""
    # Resolve file paths
    mztab_file = test_data[0]
    msstats_file = test_data[1]
    sdrf_file = test_data[2]

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


def test_extract_psm_msg():
    """Test extracting PSM messages."""
    # Resolve file paths
    mztab_file = test_data[0]
    msstats_file = test_data[1]
    sdrf_file = test_data[2]

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


def test_generate_feature():
    """Test generating features."""
    # Resolve file paths
    mztab_file = test_data[0]
    msstats_file = test_data[1]
    sdrf_file = test_data[2]

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

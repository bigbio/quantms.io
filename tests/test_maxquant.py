import pytest
from quantmsio.core.feature import Feature
from quantmsio.core.maxquant import MaxQuant
from quantmsio.core.sdrf import SDRFHandler
from quantmsio.core.psm import Psm
from .common import datafile


test_data = [
    (
        "maxquant/msms.txt",
        "maxquant/evidence.txt",
        "maxquant/sdrf.tsv",
    ),
]


@pytest.mark.parametrize("msms_path,evidence_path,sdrf_path", test_data)
def test_transform_psm(msms_path, evidence_path, sdrf_path):
    """Test transforming MaxQuant PSM data."""
    # Resolve file path
    msms_file = datafile(msms_path)

    # Initialize MaxQuant
    maxquant = MaxQuant()

    # Process PSM data
    count = 0
    for df in maxquant.iter_batch(msms_file, "psm", chunksize=500000):
        # Transform PSM data
        maxquant.transform_psm(df)
        Psm.convert_to_parquet_format(df)
        psm_parquet = Psm.transform_parquet(df)

        # Add assertions to verify the result
        assert psm_parquet is not None
        assert len(psm_parquet) > 0
        assert "peptidoform" in psm_parquet.column_names
        count += 1

    # Ensure we got at least one result
    assert count > 0


@pytest.mark.parametrize("msms_path,evidence_path,sdrf_path", test_data)
def test_transform_feature(msms_path, evidence_path, sdrf_path):
    """Test transforming MaxQuant feature data."""
    # Resolve file paths
    evidence_file = datafile(evidence_path)
    sdrf_file = datafile(sdrf_path)

    # Initialize SDRF handler and MaxQuant
    sdrf = SDRFHandler(sdrf_file)
    maxquant = MaxQuant()
    maxquant.experiment_type = sdrf.get_experiment_type_from_sdrf()
    maxquant._sample_map = sdrf.get_sample_map_run()

    # Process feature data
    count = 0
    for df in maxquant.iter_batch(evidence_file, chunksize=500000):
        # Transform feature data
        maxquant.transform_feature(df)
        Feature.convert_to_parquet_format(df)
        feature = Feature.transform_feature(df)

        # Add assertions to verify the result
        assert feature is not None
        assert len(feature) > 0
        assert "peptidoform" in feature.column_names
        count += 1

    # Ensure we got at least one result
    assert count > 0


@pytest.mark.parametrize("msms_path,evidence_path,sdrf_path", test_data)
def test_transform_features(msms_path, evidence_path, sdrf_path):
    """Test transforming MaxQuant features with slicing."""
    # Resolve file paths
    evidence_file = datafile(evidence_path)
    sdrf_file = datafile(sdrf_path)

    # Initialize SDRF handler and MaxQuant
    sdrf = SDRFHandler(sdrf_file)
    maxquant = MaxQuant()
    maxquant.experiment_type = sdrf.get_experiment_type_from_sdrf()
    maxquant._sample_map = sdrf.get_sample_map_run()

    # Process feature data
    for report in maxquant.iter_batch(evidence_file, chunksize=500000):
        # Transform feature data
        maxquant.transform_feature(report)
        Feature.convert_to_parquet_format(report)

        # Test slicing
        slice_count = 0
        for _, df in Feature.slice(report, ["reference_file_name", "precursor_charge"]):
            feature = Feature.transform_feature(df)

            # Add assertions to verify the result
            assert feature is not None
            assert len(feature) > 0
            assert "peptidoform" in feature.column_names
            slice_count += 1

        # Ensure we got at least one slice
        assert slice_count > 0

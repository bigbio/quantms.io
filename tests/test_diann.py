import pytest
from .common import datafile
from quantmsio.core.diann import DiaNNConvert
from quantmsio.core.feature import Feature


test_data = [
    (
        "DIANN/diann_report.tsv",
        "DIANN/PXD019909-DIA.sdrf.tsv",
        "DIANN/mzml",
    ),
]


@pytest.mark.parametrize("report_path,sdrf_path,mzml_path", test_data)
def test_transform_feature(report_path, sdrf_path, mzml_path):
    """Test transforming DIANN report to feature."""
    report_file = datafile(report_path)
    sdrf_file = datafile(sdrf_path)
    mzml = datafile(mzml_path)

    # Initialize DiaNNConvert
    diann = DiaNNConvert(report_file, sdrf_file)

    # Process reports
    for report in diann.main_report_df(0.05, mzml, 2):
        diann.add_additional_msg(report)
        Feature.convert_to_parquet_format(report)
        feature = Feature.transform_feature(report)

        # Add assertions to verify the result
        assert feature is not None
        assert len(feature) > 0
        assert "peptidoform" in feature.column_names


@pytest.mark.parametrize("report_path,sdrf_path,mzml_path", test_data)
def test_transform_features(report_path, sdrf_path, mzml_path):
    """Test transforming DIANN report to features with slicing."""
    report_file = datafile(report_path)
    sdrf_file = datafile(sdrf_path)
    mzml = datafile(mzml_path)

    # Initialize DiaNNConvert
    diann = DiaNNConvert(report_file, sdrf_file)

    # Process reports
    for report in diann.main_report_df(0.05, mzml, 2):
        diann.add_additional_msg(report)
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

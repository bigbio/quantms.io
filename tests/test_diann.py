from pathlib import Path

from quantmsio.core.diann import DiaNNConvert
from quantmsio.core.feature import Feature

TEST_DATA_ROOT = Path(__file__).parent / "examples"

TEST_DATA = (
    TEST_DATA_ROOT / "DIANN/diann_report.tsv",
    TEST_DATA_ROOT / "DIANN/PXD019909-DIA.sdrf.tsv",
    TEST_DATA_ROOT / "DIANN/mzml",
)


def test_transform_feature():
    report_file = TEST_DATA[0]
    sdrf_file = TEST_DATA[1]
    mzml = TEST_DATA[2]
    D = DiaNNConvert(report_file, sdrf_file)
    for report in D.main_report_df(0.05, mzml, 2):
        D.add_additional_msg(report)
        Feature.convert_to_parquet_format(report)
        Feature.transform_feature(report)


def test_transform_features():
    report_file = TEST_DATA[0]
    sdrf_file = TEST_DATA[1]
    mzml = TEST_DATA[2]
    D = DiaNNConvert(report_file, sdrf_file)
    for report in D.main_report_df(0.05, mzml, 2):
        D.add_additional_msg(report)
        Feature.convert_to_parquet_format(report)
        for _, df in Feature.slice(report, ["reference_file_name", "precursor_charge"]):
            Feature.transform_feature(df)
"""SDRF run metadata written to run.parquet."""

import pyarrow.parquet as pq
import pytest

from qpx.converters.sdrf import SdrfConverter

OXIDATION = {
    "accession": "UNIMOD:35",
    "name": "Oxidation",
    "fixed": False,
    "position": "Anywhere",
    "target_amino_acid": "M",
}


def _convert_single_run(tmp_path, modification_values):
    """Convert a one-run SDRF with one modification column per value and return the run row."""
    header = "\t".join(
        ["source name", "comment[data file]", "comment[label]"] + ["comment[modification parameters]"] * len(modification_values)
    )
    row = "\t".join(["S1", "run_01.raw", "label free sample", *modification_values])
    sdrf = tmp_path / "input.sdrf.tsv"
    sdrf.write_text(f"{header}\n{row}\n", encoding="utf-8")
    run_path = tmp_path / "run.parquet"
    with SdrfConverter() as converter:
        converter.convert(str(sdrf), run_output=str(run_path))
    return pq.read_table(run_path).to_pylist()[0]


@pytest.mark.parametrize("sentinel", ["not available", "not applicable"])
def test_sentinel_modification_columns_are_skipped(tmp_path, sentinel):
    """Sentinel modification columns must not add empty modification entries."""
    run = _convert_single_run(tmp_path, ["NT=Oxidation;MT=Variable;TA=M;AC=UNIMOD:35;PP=Anywhere", sentinel, sentinel])

    assert run["modification_parameters"] == [OXIDATION]


def test_only_sentinel_modification_columns_leave_none(tmp_path):
    """A run whose modification columns are all sentinels has no modification parameters."""
    run = _convert_single_run(tmp_path, ["not available", "not available"])

    assert run["modification_parameters"] is None

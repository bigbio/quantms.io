"""SDRF run metadata written to run.parquet."""

import pyarrow.parquet as pq
import pytest

from qpx.converters.sdrf import SdrfConverter
from qpx.core.sdrf import parse_sdrf_kv

OXIDATION = {
    "accession": "UNIMOD:35",
    "name": "Oxidation",
    "fixed": False,
    "position": "Anywhere",
    "target_amino_acid": "M",
}


def _convert_sdrf(tmp_path, columns, rows):
    """Convert an SDRF with the given header and rows; return its run rows and run ontology entries."""
    sdrf = tmp_path / "input.sdrf.tsv"
    lines = ["\t".join(columns), *("\t".join(row) for row in rows)]
    sdrf.write_text("\n".join(lines) + "\n", encoding="utf-8")
    run_path = tmp_path / "run.parquet"
    with SdrfConverter() as converter:
        converter.convert(str(sdrf), run_output=str(run_path))
        ontology_entries = converter.run_ontology_entries()
    return pq.read_table(run_path).to_pylist(), ontology_entries


def _convert_single_run(tmp_path, modification_values):
    """Convert a one-run SDRF with one modification column per value and return the run row."""
    modification_columns = ["comment[modification parameters]"] * len(modification_values)
    columns = ["source name", "comment[data file]", "comment[label]", *modification_columns]
    runs, _ = _convert_sdrf(tmp_path, columns, [["S1", "run_01.raw", "label free sample", *modification_values]])
    return runs[0]


@pytest.mark.parametrize("sentinel", ["not available", "not applicable"])
def test_sentinel_modification_columns_are_skipped(tmp_path, sentinel):
    """Sentinel modification columns must not add empty modification entries."""
    run = _convert_single_run(tmp_path, ["NT=Oxidation;MT=Variable;TA=M;AC=UNIMOD:35;PP=Anywhere", sentinel, sentinel])

    assert run["modification_parameters"] == [OXIDATION]


def test_only_sentinel_modification_columns_leave_none(tmp_path):
    """A run whose modification columns are all sentinels has no modification parameters."""
    run = _convert_single_run(tmp_path, ["not available", "not available"])

    assert run["modification_parameters"] is None


@pytest.mark.parametrize(
    ("sdrf_value", "expected"),
    [
        ("NT=Oxidation;MT=Variable;TA=M;AC=UNIMOD:35", {"NT": "Oxidation", "MT": "Variable", "TA": "M", "AC": "UNIMOD:35"}),
        (" NT = Oxidation ; AC = UNIMOD:35 ", {"NT": "Oxidation", "AC": "UNIMOD:35"}),
        ("not available", {}),
    ],
)
def test_parse_sdrf_kv_strips_pairs_and_skips_values_without_equals(sdrf_value, expected):
    """Complex SDRF values split into stripped key/value pairs, and parts without ``=`` are ignored."""
    assert parse_sdrf_kv(sdrf_value) == expected


@pytest.fixture(name="local_ms")
def _local_ms(local_ms_ontology):
    """Serve a small local PSI-MS ontology from ``PublicOntology``."""
    return local_ms_ontology(
        {
            "MS:1000133": "collision-induced dissociation",
            "MS:1000422": "beam-type collision-induced dissociation",
            "MS:1002732": "Orbitrap Fusion Lumos",
        }
    )


def _run_ontology_by_field(tmp_path, instrument, dissociation):
    """Convert a one-run SDRF and return its run ontology entries keyed by ``field_name``."""
    _, entries = _convert_sdrf(
        tmp_path,
        ["source name", "comment[data file]", "comment[label]", "comment[instrument]", "comment[dissociation method]"],
        [["S1", "run_01.raw", "label free sample", instrument, dissociation]],
    )
    return {entry["field_name"]: entry for entry in entries}


@pytest.mark.parametrize(
    ("dissociation", "field_name", "accession", "ontology_name"),
    [
        ("NT=CID;AC=MS:1000133", "CID", "MS:1000133", "collision-induced dissociation"),
        ("NT=HCD;AC=MS:1000422", "HCD", "MS:1000422", "beam-type collision-induced dissociation"),
        ("collision-induced dissociation", "collision-induced dissociation", "MS:1000133", "collision-induced dissociation"),
        ("NT=HCD;AC=MS:1000133", "HCD", "MS:1000133", "collision-induced dissociation"),
        (
            "NT=collision-induced dissociation;AC=MS:9999999",
            "collision-induced dissociation",
            "MS:1000133",
            "collision-induced dissociation",
        ),
    ],
)
def test_dissociation_resolves_by_accession_before_name(local_ms, tmp_path, dissociation, field_name, accession, ontology_name):
    """The SDRF accession wins over the name, and the name is the fallback when no accession exists."""
    entries = _run_ontology_by_field(tmp_path, instrument="not applicable", dissociation=dissociation)

    assert entries[field_name]["ontology_accession"] == accession
    assert entries[field_name]["ontology_name"] == ontology_name
    assert entries[field_name]["view"] == "run"


def test_shared_name_resolves_to_the_sorted_first_accession(local_ms, tmp_path):
    """One name with two resolvable accessions yields one entry, independent of set iteration order."""
    _, entries = _convert_sdrf(
        tmp_path,
        ["source name", "comment[data file]", "comment[label]", "comment[dissociation method]"],
        [
            ["S1", "run_01.raw", "label free sample", "NT=HCD;AC=MS:1000422"],
            ["S2", "run_02.raw", "label free sample", "NT=HCD;AC=MS:1000133"],
        ],
    )
    hcd_entries = [entry for entry in entries if entry["field_name"] == "HCD"]

    assert [entry["ontology_accession"] for entry in hcd_entries] == ["MS:1000133"]


def test_instrument_resolves_by_accession(local_ms, tmp_path):
    """An instrument ``AC=`` accession resolves its PSI-MS term when the name alone does not match a term."""
    entries = _run_ontology_by_field(
        tmp_path,
        instrument="NT=Lumos;AC=MS:1002732",
        dissociation="not applicable",
    )

    assert entries["Lumos"]["ontology_accession"] == "MS:1002732"
    assert entries["Lumos"]["ontology_name"] == "Orbitrap Fusion Lumos"

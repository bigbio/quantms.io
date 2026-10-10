"""SDRF run metadata written to run.parquet."""

import pyarrow.parquet as pq
import pytest

from qpx.converters.sdrf import SdrfConverter
from qpx.core.ontology import PublicOntology

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


@pytest.fixture(name="local_ms")
def _local_ms(tmp_path, monkeypatch):
    """Serve a small local PSI-MS ontology from ``PublicOntology``."""
    path = tmp_path / "psi-ms.obo"
    path.write_text(
        '\n[Term]\nid: MS:1000133\nname: collision-induced dissociation\ndef: "Collision-induced dissociation." []\n'
        '\n[Term]\nid: MS:1000422\nname: beam-type collision-induced dissociation\ndef: "Beam-type CID." []\n'
        '\n[Term]\nid: MS:1002732\nname: Orbitrap Fusion Lumos\ndef: "Orbitrap Fusion Lumos." []\n'
    )
    with PublicOntology.from_obo(path) as ontology:
        monkeypatch.setattr("qpx.core.ontology.PublicOntology", lambda *_args, **_kwargs: ontology)
        yield ontology


def _run_ontology_by_field(tmp_path, instrument, dissociation):
    """Convert a one-run SDRF and return its run ontology entries keyed by ``field_name``."""
    sdrf = tmp_path / "input.sdrf.tsv"
    sdrf.write_text(
        "source name\tcomment[data file]\tcomment[label]\tcomment[instrument]\tcomment[dissociation method]\n"
        f"S1\trun_01.raw\tlabel free sample\t{instrument}\t{dissociation}\n",
        encoding="utf-8",
    )
    with SdrfConverter() as converter:
        converter.convert(str(sdrf), run_output=str(tmp_path / "run.parquet"))
        return {entry["field_name"]: entry for entry in converter.run_ontology_entries()}


@pytest.mark.parametrize(
    ("dissociation", "field_name", "accession", "ontology_name"),
    [
        ("NT=CID;AC=MS:1000133", "CID", "MS:1000133", "collision-induced dissociation"),
        ("NT=HCD;AC=MS:1000422", "HCD", "MS:1000422", "beam-type collision-induced dissociation"),
        ("collision-induced dissociation", "collision-induced dissociation", "MS:1000133", "collision-induced dissociation"),
        ("NT=HCD;AC=MS:1000133", "HCD", "MS:1000133", "collision-induced dissociation"),
    ],
)
def test_dissociation_resolves_by_accession_before_name(local_ms, tmp_path, dissociation, field_name, accession, ontology_name):
    """The SDRF accession wins over the name, and the name is the fallback when no accession exists."""
    entries = _run_ontology_by_field(tmp_path, instrument="not applicable", dissociation=dissociation)

    assert entries[field_name]["ontology_accession"] == accession
    assert entries[field_name]["ontology_name"] == ontology_name
    assert entries[field_name]["view"] == "run"


def test_instrument_resolves_by_accession(local_ms, tmp_path):
    """An instrument ``AC=`` accession resolves its PSI-MS term."""
    entries = _run_ontology_by_field(
        tmp_path,
        instrument="NT=Orbitrap Fusion Lumos;AC=MS:1002732",
        dissociation="not applicable",
    )

    assert entries["Orbitrap Fusion Lumos"]["ontology_accession"] == "MS:1002732"

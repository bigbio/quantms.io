"""Modification candidates (psm/feature) and sequence edits (pepmap/proteoform)."""

import pyarrow.parquet as pq
import pytest

from tests.conftest import make_psm_record


def candidate(position, aa, name, probability, assigned, accession=None, mass=None):
    return {
        "position": position,
        "amino_acid": aa,
        "name": name,
        "accession": accession,
        "mass_delta": mass,
        "probability": probability,
        "assigned": assigned,
        "scores": None,
    }


# PEPS[Phospho]TK[Acetyl]R: phospho on S4 or T5; acetyl or trimethyl on K6.
CANDIDATES = [
    candidate(4, "S", "Phospho", 0.95, True, "UNIMOD:21", 79.966331),
    candidate(5, "T", "Phospho", 0.05, False, "UNIMOD:21", 79.966331),
    candidate(6, "K", "Acetyl", 0.71, True, "UNIMOD:1", 42.010565),
    candidate(6, "K", "Trimethyl", 0.29, False, "UNIMOD:37", 42.04695),
]


def _psm(**overrides):
    record = make_psm_record(sequence="PEPSTKR", peptidoform="PEPS[Phospho]TK[Acetyl]R")
    record["modification_candidates"] = CANDIDATES
    record.update(overrides)
    return record


def _validate(path, structure):
    import qpx

    ds = qpx.open_dataset(str(path))
    try:
        return ds.validate([structure], strict=True)[structure]
    finally:
        ds.close()


def test_candidates_round_trip_and_validate(tmp_path):
    from qpx.writers import PsmWriter

    with PsmWriter(tmp_path / "exp.psm.parquet") as w:
        w.write_batch([_psm(), make_psm_record(scan=[1002])])
    table = pq.read_table(tmp_path / "exp.psm.parquet")
    rows = table.column("modification_candidates").to_pylist()
    assert [(c["position"], c["name"], c["assigned"]) for c in rows[0]] == [
        (4, "Phospho", True),
        (5, "Phospho", False),
        (6, "Acetyl", True),
        (6, "Trimethyl", False),
    ]
    assert rows[1] is None, "a PSM without alternatives leaves the column null"
    assert _validate(tmp_path, "psm").is_valid


def test_a_file_without_the_column_stays_valid(tmp_path):
    """The column is optional: files written before it existed still validate."""
    from qpx.writers import PsmWriter

    with PsmWriter(tmp_path / "exp.psm.parquet") as w:
        w.write_batch([make_psm_record()])
    table = pq.read_table(tmp_path / "exp.psm.parquet")
    old = table.drop_columns(["modification_candidates"])
    pq.write_table(old.replace_schema_metadata(table.schema.metadata), tmp_path / "exp.psm.parquet")
    assert _validate(tmp_path, "psm").is_valid


@pytest.mark.parametrize(
    ("candidates", "check"),
    [
        ([candidate(6, "K", "Acetyl", 1.5, True)], "value_out_of_range"),
        (
            [candidate(6, "K", "Acetyl", 0.8, True), candidate(6, "K", "Trimethyl", 0.4, False)],
            "candidate_probability_sum",
        ),
    ],
)
def test_bad_candidate_probabilities_fail_strict_validation(tmp_path, candidates, check):
    from qpx.writers import PsmWriter

    with PsmWriter(tmp_path / "exp.psm.parquet") as w:
        w.write_batch([_psm(modification_candidates=candidates)])
    result = _validate(tmp_path, "psm")
    assert any(i.check == check for i in result.errors)


def edit(start, end, ref, alt, name=None):
    return {
        "protein_accession": "P04637",
        "start": start,
        "end": end,
        "reference": ref,
        "alternative": alt,
        "name": name,
        "accession": None,
    }


VALID_EDITS = [
    edit(175, 175, "R", "H", "p.Arg175His"),  # substitution
    edit(10, 11, "KE", "QK"),  # two adjacent residues replaced
    edit(30, 30, "G", ""),  # deletion
    edit(50, 49, "", "AA"),  # insertion before residue 50
]


def _pepmap(path, edits):
    from qpx.writers import PepMapWriter

    with PepMapWriter(path / "exp.pepmap.parquet") as w:
        w.write_batch(
            [
                {
                    "sequence": "VVRHCPHHER",
                    "peptidoform": "VVRHCPHHER",
                    "pg_accessions": [{"accession": "P04637", "start": 172, "end": 181, "pre": "K", "post": "C"}],
                    "sequence_edits": edits,
                    "is_unique": True,
                }
            ]
        )


def test_sequence_edits_cover_substitutions_spans_deletions_and_insertions(tmp_path):
    _pepmap(tmp_path, VALID_EDITS)
    stored = pq.read_table(tmp_path / "exp.pepmap.parquet").column("sequence_edits")[0].as_py()
    assert [(e["start"], e["end"], e["reference"], e["alternative"]) for e in stored] == [
        (175, 175, "R", "H"),
        (10, 11, "KE", "QK"),
        (30, 30, "G", ""),
        (50, 49, "", "AA"),
    ]
    assert _validate(tmp_path, "pepmap").is_valid


@pytest.mark.parametrize(
    "bad",
    [
        edit(10, 12, "KE", "QK"),  # reference shorter than its span
        edit(0, 0, "", "M"),  # start before the protein
        edit(20, 20, "r", "h"),  # not one-letter codes
        edit(20, 19, "", ""),  # changes nothing
    ],
)
def test_malformed_sequence_edits_fail_strict_validation(tmp_path, bad):
    _pepmap(tmp_path, [bad])
    result = _validate(tmp_path, "pepmap")
    assert any(i.check == "malformed_sequence_edit" for i in result.errors)


def test_a_variant_proteoform_carries_its_edits(tmp_path):
    from qpx.writers import ProteoformWriter

    record = {
        "family": "P04637",
        "accession": "P04637_VAR_R175H",
        "kind": "variant",
        "description": "p.Arg175His",
        "edits": [{"start": 175, "end": 175, "reference": "R", "alternative": "H", "name": "p.Arg175His", "accession": None}],
        "grouped_runs": ["run_01"],
        "is_decoy": False,
    }
    with ProteoformWriter(tmp_path / "exp.proteoform.parquet") as w:
        w.write_batch([record])
    stored = pq.read_table(tmp_path / "exp.proteoform.parquet").column("edits")[0].as_py()
    assert stored[0]["reference"] == "R" and stored[0]["alternative"] == "H"
    assert _validate(tmp_path, "proteoform").is_valid


def test_maxquant_phospho_probabilities_keep_every_candidate_residue():
    """S 0.95 is placed; T 0.03 and Y 0.02 were dropped before and are now candidates."""
    from qpx.converters.maxquant.constants import parse_phospho_probabilities, phospho_candidates
    from qpx.converters.ptm import from_proforma

    scores = parse_phospho_probabilities("_PEPS(0.95)T(0.03)Y(0.02)K_")
    _, modifications = from_proforma("PEPS[UNIMOD:21]TYK", "PEPSTYK", site_scores=scores)
    candidates = phospho_candidates(scores, "PEPSTYK", modifications)
    assert [(c["position"], c["amino_acid"], c["probability"], c["assigned"]) for c in candidates] == [
        (4, "S", 0.95, True),
        (5, "T", 0.03, False),
        (6, "Y", 0.02, False),
    ]
    assert phospho_candidates(None, "PEPSTYK", modifications) is None

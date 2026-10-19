"""Site and proteoform view tests."""

import pyarrow.parquet as pq
import pytest


def make_site_record(**overrides):
    defaults = {
        "protein_accession": "P04637",
        "family": "P04637",
        "position": 120,
        "amino_acid": "K",
        "modification_name": "Acetyl",
        "modification_accession": "UNIMOD:1",
        "mass_delta": 42.010565,
        "grouped_runs": ["run_01"],
        "pg_id": None,
        "prior_evidence": "curated",
        "site_qvalue": 0.004,
        "occupancy": 0.97,
        "identity_probability": 0.71,
        "alternatives": [
            {"name": "Trimethyl", "accession": "UNIMOD:37", "mass_delta": 42.04695, "identity_probability": 0.29},
        ],
        "psm_count": 12,
        "best_psm_id": None,
        "is_decoy": False,
    }
    defaults.update(overrides)
    return defaults


def make_proteoform_record(**overrides):
    defaults = {
        "family": "P04637",
        "accession": "P04637-VAR_005932",
        "kind": "variant",
        "description": "p.Arg72Pro",
        "pg_id": None,
        "grouped_runs": ["run_01"],
        "source": "UniProt",
        "prior_evidence": "curated",
        "specific_peptides": 2,
        "shared_peptides": 9,
        "psm_count": 5,
        "evidence": "specific",
        "qvalue": 0.01,
        "best_psm_id": None,
        "is_decoy": False,
    }
    defaults.update(overrides)
    return defaults


def test_site_schema_round_trip_and_identity(tmp_path):
    """Two modifications claimed at one residue are two rows with distinct ids."""
    from qpx.core.data import SiteSchema
    from qpx.writers import SiteWriter

    names = [f.name for f in SiteSchema.get_arrow_schema()]
    for field in ("site_id", "protein_accession", "position", "mass_delta", "occupancy", "identity_probability"):
        assert field in names

    records = [
        make_site_record(),
        make_site_record(
            modification_name="Trimethyl",
            modification_accession="UNIMOD:37",
            mass_delta=42.04695,
            identity_probability=0.29,
            alternatives=[{"name": "Acetyl", "accession": "UNIMOD:1", "mass_delta": 42.010565, "identity_probability": 0.71}],
        ),
        # A site found without a prior and without a site-level estimate.
        make_site_record(position=382, prior_evidence=None, site_qvalue=None, alternatives=None),
    ]
    out = tmp_path / "exp.site.parquet"
    with SiteWriter(out) as w:
        w.write_batch(records)
    table = pq.read_table(out)
    assert table.num_rows == 3
    ids = table.column("site_id").to_pylist()
    assert len(set(ids)) == 3, "mass_delta is part of the identity"
    alt = table.column("alternatives")[0].as_py()
    assert alt[0]["name"] == "Trimethyl"
    assert alt[0]["identity_probability"] == pytest.approx(0.29)
    assert table.column("alternatives")[2].as_py() is None


def test_proteoform_rows_per_label(tmp_path):
    """Quantities unroll into one row per label, as in pg; ID-only members stay one row."""
    from qpx.writers import ProteoformWriter

    records = [
        make_proteoform_record(
            intensities=[{"label": "TMT126", "intensity": 10.0}, {"label": "TMT127N", "intensity": 12.5}],
        ),
        make_proteoform_record(accession="P04637", kind="canonical", description=None, specific_peptides=20),
    ]
    out = tmp_path / "exp.proteoform.parquet"
    with ProteoformWriter(out) as w:
        w.write_batch(records)
    table = pq.read_table(out)
    assert table.num_rows == 3
    labels = sorted((r or "") for r in table.column("label").to_pylist())
    assert labels == ["", "TMT126", "TMT127N"]
    assert len(set(table.column("proteoform_id").to_pylist())) == 3


def test_dataset_discovers_site_and_proteoform(tmp_path):
    """Dataset finds both views; their filters work."""
    import qpx
    from qpx.writers import ProteoformWriter, SiteWriter

    with SiteWriter(tmp_path / "exp.site.parquet") as w:
        w.write_batch(
            [
                make_site_record(),
                make_site_record(position=382, site_qvalue=0.2, occupancy=0.4),
                make_site_record(protein_accession="P38398", family="P38398", modification_name="Phospho", is_decoy=True),
            ]
        )
    with ProteoformWriter(tmp_path / "exp.proteoform.parquet") as w:
        w.write_batch(
            [
                make_proteoform_record(),
                make_proteoform_record(accession="P04637-2", kind="isoform", description="isoform 2", qvalue=0.2),
                make_proteoform_record(accession="P04637", kind="canonical", description=None),
            ]
        )

    ds = qpx.open_dataset(str(tmp_path))
    assert ds.site.count() == 3
    assert ds.site.by_protein("P04637").count() == 2
    assert ds.site.by_modification("Phospho").count() == 1
    assert ds.site.targets_only().confident(0.05, min_occupancy=0.5, min_identity=0.7).count() == 1
    assert ds.proteoform.by_family("P04637").count() == 3
    assert ds.proteoform.by_kind("variant").count() == 1
    assert ds.proteoform.confident(0.05).count() == 2
    with pytest.raises(ValueError):
        ds.proteoform.by_kind("gene")
    ds.close()

    assert qpx.read_site(str(tmp_path / "exp.site.parquet")).count() == 3
    assert qpx.read_proteoform(str(tmp_path / "exp.proteoform.parquet")).count() == 3

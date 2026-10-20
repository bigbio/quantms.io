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
        "modified_probability": 0.97,
        "evidence_fraction": 0.71,
        "alternatives": [
            {"name": "Trimethyl", "accession": "UNIMOD:37", "mass_delta": 42.04695, "evidence_fraction": 0.29},
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
        "proteoform_qvalue": 0.01,
        "best_psm_id": None,
        "is_decoy": False,
    }
    defaults.update(overrides)
    return defaults


def _validate(tmp_path, structure):
    import qpx

    ds = qpx.open_dataset(str(tmp_path))
    try:
        return ds.validate([structure], strict=True)[structure]
    finally:
        ds.close()


def test_site_schema_round_trip_and_identity(tmp_path):
    """Two modifications claimed at one residue are two rows with distinct ids."""
    from qpx.core.data import SiteSchema
    from qpx.writers import SiteWriter

    names = [f.name for f in SiteSchema.get_arrow_schema()]
    for field in ("site_id", "protein_accession", "position", "mass_delta", "modified_probability", "evidence_fraction"):
        assert field in names
    for renamed in ("occupancy", "identity_probability"):
        assert renamed not in names

    records = [
        make_site_record(),
        make_site_record(
            modification_name="Trimethyl",
            modification_accession="UNIMOD:37",
            mass_delta=42.04695,
            evidence_fraction=0.29,
            alternatives=[{"name": "Acetyl", "accession": "UNIMOD:1", "mass_delta": 42.010565, "evidence_fraction": 0.71}],
        ),
        # A site found without a prior and without a site-level estimate.
        make_site_record(position=382, prior_evidence=None, site_qvalue=None, alternatives=None, evidence_fraction=1.0),
    ]
    out = tmp_path / "exp.site.parquet"
    with SiteWriter(out) as w:
        w.write_batch(records)
    table = pq.read_table(out)
    assert table.num_rows == 3
    assert len(set(table.column("site_id").to_pylist())) == 3
    alt = table.column("alternatives")[0].as_py()
    assert alt[0]["name"] == "Trimethyl"
    assert alt[0]["evidence_fraction"] == pytest.approx(0.29)
    assert table.column("alternatives")[2].as_py() is None


def test_equal_mass_modifications_get_distinct_site_ids(tmp_path):
    """Dimethyl and ethyl (both C2H4) at one residue are two claims, two ids; strict validation passes."""
    from qpx.writers import SiteWriter

    with SiteWriter(tmp_path / "exp.site.parquet") as w:
        w.write_batch(
            [
                make_site_record(
                    modification_name="Dimethyl",
                    modification_accession="UNIMOD:36",
                    mass_delta=28.0313,
                    alternatives=None,
                    evidence_fraction=1.0,
                ),
                make_site_record(
                    modification_name="Ethyl",
                    modification_accession="UNIMOD:280",
                    mass_delta=28.0313,
                    alternatives=None,
                    evidence_fraction=1.0,
                ),
            ]
        )
    ids = pq.read_table(tmp_path / "exp.site.parquet").column("site_id").to_pylist()
    assert ids[0] != ids[1]
    assert _validate(tmp_path, "site").is_valid


def test_protein_n_terminus_and_residue_1_are_distinct_sites(tmp_path):
    """N-terminal acetylation is position 0 with no residue; a Lys1 side chain is position 1."""
    from qpx.writers import SiteWriter

    with SiteWriter(tmp_path / "exp.site.parquet") as w:
        w.write_batch(
            [
                make_site_record(position=0, amino_acid=None, alternatives=None, evidence_fraction=1.0),
                make_site_record(position=1, amino_acid="K", alternatives=None, evidence_fraction=1.0),
            ]
        )
    table = pq.read_table(tmp_path / "exp.site.parquet")
    assert table.column("amino_acid").to_pylist() == [None, "K"]
    assert len(set(table.column("site_id").to_pylist())) == 2
    assert _validate(tmp_path, "site").is_valid


@pytest.mark.parametrize(
    "overrides",
    [
        {"site_qvalue": -0.5},
        {"modified_probability": 1.5},
        {"evidence_fraction": 1.2, "alternatives": None},
        {"evidence_fraction": float("nan"), "alternatives": None},
        {
            "evidence_fraction": 0.5,
            "alternatives": [{"name": "Trimethyl", "accession": None, "mass_delta": 42.04695, "evidence_fraction": 1.5}],
        },
        # Without alternatives the row holds the residue's whole evidence.
        {"evidence_fraction": 0.71, "alternatives": None},
        # Alternatives listed smallest first.
        {
            "evidence_fraction": 0.5,
            "alternatives": [
                {"name": "Methyl", "accession": None, "mass_delta": 14.01565, "evidence_fraction": 0.1},
                {"name": "Trimethyl", "accession": None, "mass_delta": 42.04695, "evidence_fraction": 0.4},
            ],
        },
        # In range, but the residue's fractions do not sum to 1.
        {
            "evidence_fraction": 0.71,
            "alternatives": [{"name": "Trimethyl", "accession": None, "mass_delta": 42.04695, "evidence_fraction": 0.5}],
        },
    ],
)
def test_site_values_outside_their_range_fail_strict_validation(tmp_path, overrides):
    from qpx.writers import SiteWriter

    with SiteWriter(tmp_path / "exp.site.parquet") as w:
        w.write_batch([make_site_record(**overrides)])
    result = _validate(tmp_path, "site")
    assert not result.is_valid
    assert {i.check for i in result.errors} & {"value_out_of_range", "evidence_fraction_sum", "alternatives_order"}


def test_proteoform_qvalue_outside_its_range_fails_strict_validation(tmp_path):
    from qpx.writers import ProteoformWriter

    with ProteoformWriter(tmp_path / "exp.proteoform.parquet") as w:
        w.write_batch([make_proteoform_record(proteoform_qvalue=2.0)])
    result = _validate(tmp_path, "proteoform")
    assert any(i.check == "value_out_of_range" for i in result.errors)


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


def _two_channel_pg(tmp_path):
    """One protein group quantified in TMT126 and TMT127N; returns {label: pg_id}."""
    from qpx.writers import PgWriter
    from tests.conftest import make_pg_record

    with PgWriter(tmp_path / "exp.pg.parquet") as w:
        w.write_batch(
            [
                make_pg_record(
                    anchor_protein="P04637",
                    intensities=[{"label": "TMT126", "intensity": 1.0}, {"label": "TMT127N", "intensity": 2.0}],
                )
            ]
        )
    pg = pq.read_table(tmp_path / "exp.pg.parquet")
    return dict(zip(pg.column("label").to_pylist(), pg.column("pg_id").to_pylist()))


def test_proteoform_links_each_label_to_its_own_pg_row(tmp_path):
    """pg_id per intensities entry: each label row joins the pg row of that label."""
    from qpx.writers import ProteoformWriter

    pg_ids = _two_channel_pg(tmp_path)
    with ProteoformWriter(tmp_path / "exp.proteoform.parquet") as w:
        w.write_batch(
            [
                make_proteoform_record(
                    intensities=[
                        {"label": "TMT126", "intensity": 10.0, "pg_id": pg_ids["TMT126"]},
                        {"label": "TMT127N", "intensity": 12.5, "pg_id": pg_ids["TMT127N"]},
                    ]
                )
            ]
        )
    table = pq.read_table(tmp_path / "exp.proteoform.parquet")
    links = dict(zip(table.column("label").to_pylist(), table.column("pg_id").to_pylist()))
    assert links == pg_ids
    assert _validate(tmp_path, "proteoform").is_valid


def test_a_record_level_pg_id_for_several_labels_is_refused(tmp_path):
    from qpx.writers import ProteoformWriter

    pg_ids = _two_channel_pg(tmp_path)
    with pytest.raises(ValueError, match="pg_id"):
        with ProteoformWriter(tmp_path / "exp.proteoform.parquet") as w:
            w.write_batch(
                [
                    make_proteoform_record(
                        pg_id=pg_ids["TMT126"],
                        intensities=[{"label": "TMT126", "intensity": 10.0}, {"label": "TMT127N", "intensity": 12.5}],
                    )
                ]
            )


def test_a_cross_label_pg_link_fails_strict_validation(tmp_path):
    """A flat TMT127N row that names the TMT126 pg row is reported."""
    from qpx.writers import ProteoformWriter

    pg_ids = _two_channel_pg(tmp_path)
    with ProteoformWriter(tmp_path / "exp.proteoform.parquet") as w:
        w.write_batch([make_proteoform_record(label="TMT127N", intensity=12.5, pg_id=pg_ids["TMT126"])])
    result = _validate(tmp_path, "proteoform")
    assert any(i.check == "pg_link_context" for i in result.errors)


@pytest.mark.parametrize("overrides", [{"kind": "gene"}, {"evidence": "maybe"}])
def test_proteoform_values_outside_the_vocabulary_fail_strict_validation(tmp_path, overrides):
    from qpx.writers import ProteoformWriter

    with ProteoformWriter(tmp_path / "exp.proteoform.parquet") as w:
        w.write_batch([make_proteoform_record(**overrides)])
    result = _validate(tmp_path, "proteoform")
    assert any(i.check == "value_not_in_vocabulary" for i in result.errors)


def _id_only_pg(tmp_path):
    """One identification-only protein group (null label); returns its pg_id."""
    from qpx.writers import PgWriter
    from tests.conftest import make_pg_record

    with PgWriter(tmp_path / "exp.pg.parquet") as w:
        w.write_batch([make_pg_record(anchor_protein="P04637", intensities=[])])
    pg = pq.read_table(tmp_path / "exp.pg.parquet")
    assert pg.column("label").to_pylist() == [None]
    return pg.column("pg_id")[0].as_py()


def test_a_site_links_to_the_identification_pg_row(tmp_path):
    """A site links to the null-label pg row of its runs; a labelled pg row is reported."""
    from qpx.writers import SiteWriter

    pg_id = _id_only_pg(tmp_path)
    with SiteWriter(tmp_path / "exp.site.parquet") as w:
        w.write_batch([make_site_record(pg_id=pg_id)])
    assert _validate(tmp_path, "site").is_valid

    for f in tmp_path.glob("exp.*.parquet"):
        f.unlink()
    labelled = _two_channel_pg(tmp_path)
    with SiteWriter(tmp_path / "exp.site.parquet") as w:
        w.write_batch([make_site_record(pg_id=labelled["TMT126"])])
    result = _validate(tmp_path, "site")
    assert any(i.check == "pg_link_context" for i in result.errors)


def test_best_psm_id_must_resolve_when_the_psm_view_is_present(tmp_path):
    from qpx.writers import ProteoformWriter, PsmWriter, SiteWriter
    from tests.conftest import make_psm_record

    with PsmWriter(tmp_path / "exp.psm.parquet") as w:
        w.write_batch([make_psm_record()])
    psm_id = pq.read_table(tmp_path / "exp.psm.parquet").column("psm_id")[0].as_py()
    with SiteWriter(tmp_path / "exp.site.parquet") as w:
        w.write_batch([make_site_record(best_psm_id=psm_id), make_site_record(position=200, best_psm_id=psm_id + 1)])
    with ProteoformWriter(tmp_path / "exp.proteoform.parquet") as w:
        w.write_batch([make_proteoform_record(best_psm_id=psm_id)])
    site = _validate(tmp_path, "site")
    assert [i.check for i in site.errors] == ["dangling_psm_id"]
    assert "1 site row(s)" in site.errors[0].message
    assert _validate(tmp_path, "proteoform").is_valid


def test_dataset_discovers_site_and_proteoform(tmp_path):
    """Dataset finds both views; their filters work."""
    import qpx
    from qpx.writers import ProteoformWriter, SiteWriter

    with SiteWriter(tmp_path / "exp.site.parquet") as w:
        w.write_batch(
            [
                make_site_record(),
                make_site_record(position=382, site_qvalue=0.2, modified_probability=0.4),
                make_site_record(protein_accession="P38398", family="P38398", modification_name="Phospho", is_decoy=True),
            ]
        )
    with ProteoformWriter(tmp_path / "exp.proteoform.parquet") as w:
        w.write_batch(
            [
                make_proteoform_record(),
                make_proteoform_record(accession="P04637-2", kind="isoform", description="isoform 2", proteoform_qvalue=0.2),
                # Shared-only support, no member-level estimate.
                make_proteoform_record(
                    accession="P04637", kind="canonical", description=None, evidence="shared", proteoform_qvalue=None
                ),
            ]
        )

    ds = qpx.open_dataset(str(tmp_path))
    assert ds.site.count() == 3
    assert ds.site.by_protein("P04637").count() == 2
    assert ds.site.by_modification("Phospho").count() == 1
    assert ds.site.targets_only().confident(0.05, min_modified_probability=0.5, min_evidence_fraction=0.7).count() == 1
    assert ds.site.confident(None).count() == 3
    assert ds.proteoform.by_family("P04637").count() == 3
    assert ds.proteoform.by_kind("variant").count() == 1
    assert ds.proteoform.confident(0.05).count() == 1
    with pytest.raises(ValueError):
        ds.proteoform.by_kind("gene")
    ds.close()

    assert qpx.read_site(str(tmp_path / "exp.site.parquet")).count() == 3
    assert qpx.read_proteoform(str(tmp_path / "exp.proteoform.parquet")).count() == 3

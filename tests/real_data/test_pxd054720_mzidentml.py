"""Real-data integration test — PXD054720 mzIdentML (Mascot 1.2).

Converts the compressed mzIdentML file from PRIDE project PXD054720 to QPX
PSM Parquet and validates the output against schema expectations.

Files:
    tests/examples/mzidentml/PXD054720/F001234.mzid.gz  (~1.2 MB, 9 707 spectra)
    tests/examples/mzidentml/PXD054720/F001234.mgf.gz   (~25 MB, raw spectra)

Markers:
    @pytest.mark.large_data — CI runs it with the integration tests, not the unit tests
"""

import pytest

lxml = pytest.importorskip("lxml", reason="lxml required for mzIdentML tests")

from pathlib import Path

import pyarrow.parquet as pq

from qpx.converters.mzidentml.psm_adapter import MzIdentMLPsmAdapter
from qpx.core.data import PsmSchema

EXAMPLES = Path(__file__).parent.parent / "examples" / "mzidentml" / "PXD054720"
MZID_GZ = EXAMPLES / "F001234.mzid.gz"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def psm_table(tmp_path_factory):
    """Convert the compressed mzid.gz directly to PSM Parquet once per module."""
    if not MZID_GZ.exists():
        pytest.skip(f"Test data not found: {MZID_GZ}")

    work_dir = tmp_path_factory.mktemp("pxd054720")

    # Convert directly from .mzid.gz (adapter handles gzip transparently)
    out_path = work_dir / "F001234.psm.parquet"
    with MzIdentMLPsmAdapter() as adapter:
        adapter.convert(mzid_path=str(MZID_GZ), output_path=str(out_path))

    return pq.read_table(str(out_path))


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestPXD054720Conversion:
    """Validate MzIdentML conversion on real PRIDE data."""

    def test_produces_psms(self, psm_table):
        """Conversion must produce a non-trivial number of PSMs."""
        assert psm_table.num_rows > 1000, f"Expected >1000 PSMs from 9707 spectra, got {psm_table.num_rows}"

    def test_schema_validation(self, psm_table):
        """Output must pass QPX PsmSchema validation."""
        errors = PsmSchema.validate(psm_table)
        assert not errors, f"Schema validation errors: {errors}"

    def test_required_fields_present(self, psm_table):
        """All required PSM fields must be present."""
        required = [
            "sequence",
            "peptidoform",
            "charge",
            "is_decoy",
            "calculated_mz",
            "observed_mz",
            "run_file_name",
            "scan",
            "protein_accessions",
        ]
        for field in required:
            assert field in psm_table.column_names, f"Missing required field: {field}"

    def test_run_file_name(self, psm_table):
        """All PSMs should reference the source run file."""
        runs = set(psm_table.column("run_file_name").to_pylist())
        assert len(runs) >= 1, "Expected at least one run file name"
        # The mzid references F001234.mgf → run name should be F001234
        assert "F001234" in runs or any("F001234" in r for r in runs)

    def test_sequences_are_valid(self, psm_table):
        """Peptide sequences should be non-empty and contain valid AA characters."""
        seqs = psm_table.column("sequence").to_pylist()
        assert all(len(s) > 0 for s in seqs), "Found empty sequences"
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        for seq in seqs[:100]:  # spot-check first 100
            assert set(seq).issubset(valid_aa), f"Invalid characters in sequence: {seq}"

    def test_charges_in_range(self, psm_table):
        """Charge states should be positive and reasonable."""
        charges = psm_table.column("charge").to_pylist()
        assert all(1 <= c <= 10 for c in charges), f"Unexpected charge values: {set(charges)}"

    def test_mz_values_positive(self, psm_table):
        """Both calculated and observed m/z should be positive."""
        calc_mz = psm_table.column("calculated_mz").to_pylist()
        obs_mz = psm_table.column("observed_mz").to_pylist()
        assert all(v > 0 for v in calc_mz if v is not None)
        assert all(v > 0 for v in obs_mz if v is not None)

    def test_scan_numbers(self, psm_table):
        """Scan numbers should be non-empty lists of positive integers."""
        scans = [row.as_py() for row in psm_table.column("scan")]
        assert all(len(s) > 0 for s in scans), "Found empty scan lists"
        assert all(s[0] >= 0 for s in scans), "Found negative scan numbers"

    def test_protein_accessions_populated(self, psm_table):
        """Most PSMs should have at least one protein accession."""
        prots = psm_table.column("protein_accessions").to_pylist()
        with_proteins = sum(1 for p in prots if p and len(p) > 0)
        ratio = with_proteins / len(prots)
        assert ratio > 0.9, f"Only {ratio:.1%} of PSMs have protein accessions (expected >90%)"

    def test_crosslinks_field_present(self, psm_table):
        """Cross-links column must exist (may or may not have entries)."""
        assert "cross_links" in psm_table.column_names

    def test_has_scores(self, psm_table):
        """PSMs should have search engine scores."""
        scores_col = psm_table.column("additional_scores").to_pylist()
        with_scores = sum(1 for s in scores_col if s and len(s) > 0)
        ratio = with_scores / len(scores_col)
        assert ratio > 0.5, f"Only {ratio:.1%} of PSMs have scores (expected >50%)"

    def test_unique_peptides(self, psm_table):
        """Should identify a reasonable number of unique peptides."""
        seqs = set(psm_table.column("sequence").to_pylist())
        assert len(seqs) > 100, f"Expected >100 unique peptides, got {len(seqs)}"

    def test_unique_proteins(self, psm_table):
        """Should map to a reasonable number of proteins."""
        all_prots = set()
        for p in psm_table.column("protein_accessions").to_pylist():
            if p:
                all_prots.update(p)
        assert len(all_prots) > 10, f"Expected >10 unique proteins, got {len(all_prots)}"

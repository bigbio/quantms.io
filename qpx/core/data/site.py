"""Site data structure — protein modification sites."""

from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema
from qpx.core.query import _escape_sql_string

SiteSchema = load_schema("site")


class Site(BaseStructure):
    """Protein modification-site claims with their q-value, modified probability and evidence fraction."""

    _schema_class = SiteSchema

    def by_protein(self, protein: str) -> "Site":
        """Filter sites on a protein accession."""
        return self.filter(f"protein_accession = '{_escape_sql_string(protein)}'")

    def by_modification(self, name: str) -> "Site":
        """Filter sites of one modification (by name, e.g. Phospho)."""
        return self.filter(f"modification_name = '{_escape_sql_string(name)}'")

    def confident(
        self,
        max_qvalue: float | None = 0.05,
        min_modified_probability: float | None = None,
        min_evidence_fraction: float | None = None,
    ) -> "Site":
        """Sites at a site q-value cutoff and, optionally, floors on the
        probability that the residue is modified and on this modification's
        evidence fraction.

        ``max_qvalue=None`` skips the q-value cutoff, for producers that do not
        estimate ``site_qvalue``. A row with a null value fails any floor set
        on it.
        """
        clauses = []
        if max_qvalue is not None:
            clauses.append(f"site_qvalue <= {float(max_qvalue)}")
        if min_modified_probability is not None:
            clauses.append(f"modified_probability >= {float(min_modified_probability)}")
        if min_evidence_fraction is not None:
            clauses.append(f"evidence_fraction >= {float(min_evidence_fraction)}")
        return self.filter(" AND ".join(clauses)) if clauses else self

    def targets_only(self) -> "Site":
        """Filter to target sites (exclude decoy and knockoff sites)."""
        return self.filter("is_decoy = false")

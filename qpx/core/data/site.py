"""Site data structure — protein modification sites."""

from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema
from qpx.core.query import _escape_sql_string

SiteSchema = load_schema("site")


class Site(BaseStructure):
    """Protein modification-site claims with occupancy and identity probabilities."""

    _schema_class = SiteSchema

    def by_protein(self, protein: str) -> "Site":
        """Filter sites on a protein accession."""
        return self.filter(f"protein_accession = '{_escape_sql_string(protein)}'")

    def by_modification(self, name: str) -> "Site":
        """Filter sites of one modification (by name, e.g. Phospho)."""
        return self.filter(f"modification_name = '{_escape_sql_string(name)}'")

    def confident(
        self,
        max_qvalue: float = 0.05,
        min_occupancy: float | None = None,
        min_identity: float | None = None,
    ) -> "Site":
        """Sites at a site q-value cutoff and, optionally, occupancy and identity floors."""
        clauses = [f"site_qvalue <= {float(max_qvalue)}"]
        if min_occupancy is not None:
            clauses.append(f"occupancy >= {float(min_occupancy)}")
        if min_identity is not None:
            clauses.append(f"identity_probability >= {float(min_identity)}")
        return self.filter(" AND ".join(clauses))

    def targets_only(self) -> "Site":
        """Filter to target sites (exclude decoy and knockoff sites)."""
        return self.filter("is_decoy = false")

"""PepMap data structure — peptide-to-protein mapping."""

from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema
from qpx.core.query import _escape_sql_string

PepMapSchema = load_schema("pepmap")


class PepMap(BaseStructure):
    """Deduplicated peptide-to-protein mapping."""

    _schema_class = PepMapSchema

    def by_protein(self, protein: str) -> "PepMap":
        """Filter mappings that include a given protein accession."""
        escaped = _escape_sql_string(protein)
        return self.filter(f"len(list_filter(pg_accessions, lambda x: x.accession = '{escaped}')) > 0")

    def by_peptide(self, sequence: str) -> "PepMap":
        """Filter mappings by peptide sequence."""
        return self.filter(f"sequence = '{_escape_sql_string(sequence)}'")

    def unique_peptides(self) -> "PepMap":
        """Filter to peptides that map to exactly one protein."""
        return self.filter("is_unique = true")

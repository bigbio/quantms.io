"""Proteoform data structure — family members with their own evidence."""

from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema
from qpx.core.query import _escape_sql_string

ProteoformSchema = load_schema("proteoform")

KINDS = ("canonical", "isoform", "variant", "chain")


class Proteoform(BaseStructure):
    """Canonical, isoform, variant and chain members of protein families."""

    _schema_class = ProteoformSchema

    def by_family(self, family: str) -> "Proteoform":
        """Filter members of one protein family."""
        return self.filter(f"family = '{_escape_sql_string(family)}'")

    def by_kind(self, kind: str) -> "Proteoform":
        """Filter members of one kind: canonical, isoform, variant or chain."""
        if kind not in KINDS:
            raise ValueError(f"kind must be one of {KINDS}, got {kind!r}")
        return self.filter(f"kind = '{kind}'")

    def confident(self, max_qvalue: float = 0.05) -> "Proteoform":
        """Members at a member-level q-value cutoff (rows without one are left out)."""
        return self.filter(f"proteoform_qvalue <= {float(max_qvalue)}")

    def targets_only(self) -> "Proteoform":
        """Filter to target members (exclude decoys and knockoffs)."""
        return self.filter("is_decoy = false")

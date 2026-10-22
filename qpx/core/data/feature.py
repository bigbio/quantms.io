"""Feature data structure — quantified peptide features per MS run."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from qpx.core.convert import QueryResult
from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema
from qpx.core.query import _escape_sql_string
from qpx.core.sql import sql_build

logger = logging.getLogger(__name__)

FeatureSchema = load_schema("feature")

if TYPE_CHECKING:
    from qpx.dataset import Dataset


class Feature(BaseStructure):
    """Quantified peptide features per MS run."""

    _schema_class = FeatureSchema

    def unique_proteins(self) -> list[str]:
        """Get all unique anchor proteins."""
        stmt = sql_build("SELECT DISTINCT anchor_protein FROM $src", src=self._query.source)
        result = self._engine.execute(stmt)
        return [row[0] for row in result.fetchall()]

    def by_protein(self, protein: str) -> Feature:
        """Filter features by anchor protein."""
        return self.filter(f"anchor_protein = '{_escape_sql_string(protein)}'")

    def by_run(self, run_file_name: str) -> Feature:
        """Filter features by run file."""
        return self.filter(f"run_file_name = '{_escape_sql_string(run_file_name)}'")

    def peptide_intensities(self) -> QueryResult:
        """
        Flatten intensities: one row per (peptide, label, intensity).
        Useful for downstream analysis and plotting.
        """
        ilf = self._intensity_label_field()
        stmt = sql_build(
            """SELECT sequence, peptidoform, charge, anchor_protein,
               run_file_name, i.$ilf AS label, i.intensity
        FROM $src,
             UNNEST(intensities) AS _t(i)""",
            ilf=ilf,
            src=self._query.source,
        )
        return QueryResult(self._engine.connection, stmt)

    def _intensity_label_field(self) -> str:
        """Detect whether the intensities struct uses 'label' (new) or 'channel' (old)."""
        try:
            stmt = sql_build("SELECT typeof(intensities) FROM $src LIMIT 1", src=self._query.source)
            rows = self._engine.execute(stmt).fetchall()
            if rows:
                type_str = rows[0][0].lower()
                if "channel" in type_str and "label" not in type_str:
                    return "channel"
        except Exception:
            logger.debug("Failed to detect intensities struct type", exc_info=True)
        return "label"

    def for_quantification(self, dataset: Dataset) -> QueryResult:
        """
        Return feature intensities in long format with sample metadata,
        ready for downstream quantification tools.

        Returns one row per (peptidoform, charge, protein, sample, label)
        with columns: sequence, peptidoform, charge, anchor_protein,
        sample_accession, label, intensity, run_file_name,
        biological_replicate, fraction.
        """
        ilf = self._intensity_label_field()
        # Old schema has sample_accession in intensities struct, join on that.
        # New schema has label in intensities struct, join on label.
        if ilf == "channel":
            join_cond = "i.sample_accession = rs.sample_accession"
            label_expr = "i.channel"
        else:
            join_cond = "(i.label = rs.label OR len(r.samples) = 1)"
            label_expr = "rs.label"
        stmt = sql_build(
            """SELECT f.sequence, f.peptidoform, f.charge, f.anchor_protein,
               rs.sample_accession, $label_expr AS label, i.intensity,
               f.run_file_name, rs.biological_replicate, r.fraction
        FROM $src f
        JOIN run r ON f.run_file_name = r.run_file_name,
             UNNEST(r.samples) AS _t1(rs),
             UNNEST(f.intensities) AS _t2(i)
        WHERE $join_cond""",
            ilf=ilf,
            src=self._query.source,
            join_cond=join_cond,
            label_expr=label_expr,
        )
        return QueryResult(self._engine.connection, stmt)

"""ProteoformWriter — writes proteoform.parquet files.

The proteoform view is flattened like pg: one row per label. Records may carry
the natural ``intensities`` list (``list<{label, intensity}>``); the writer
unrolls it with the same rule as :mod:`qpx.writers.pg`, and an
identification-only member becomes one row with ``label``/``intensity`` null.

``pg_id`` names one pg row, and a pg row has one label. A member quantified in
several labels therefore links each label row to that label's pg row: give
``pg_id`` on each ``intensities`` entry. A record-level ``pg_id`` is accepted
only when the record has at most one label.
"""

from __future__ import annotations

import pyarrow as pa

from qpx.core.data import ProteoformSchema
from qpx.writers.pg import PgWriter


def _split_label_links(records: list[dict]) -> list[dict]:
    """Give each label its own pg_id before the per-label explode.

    A record whose ``intensities`` entries carry ``pg_id`` is split into one
    record per entry, each with that entry's ``pg_id``. A record-level
    ``pg_id`` on a record with several labels is refused: it would link every
    label row to one label's pg row.
    """
    out: list[dict] = []
    for record in records:
        intensities = record.get("intensities") or []
        entries = [e for e in intensities if isinstance(e, dict)]
        per_entry = any("pg_id" in e for e in entries)
        if not per_entry:
            if record.get("pg_id") is not None and len(entries) > 1:
                raise ValueError(
                    "proteoform record has a record-level pg_id and "
                    f"{len(entries)} labels; pg_id names one pg row (one label): "
                    "give pg_id on each intensities entry instead"
                )
            out.append(record)
            continue
        for entry in entries:
            single = dict(record)
            single["pg_id"] = entry.get("pg_id")
            single["intensities"] = [{k: v for k, v in entry.items() if k != "pg_id"}]
            out.append(single)
    return out


class ProteoformWriter(PgWriter):
    """Write proteoform.parquet files, one row per label as in pg."""

    _schema_class = ProteoformSchema

    def _write_batch(self, records: list[dict]):
        """Resolve per-label pg links, then explode as pg does."""
        super()._write_batch(_split_label_links(records))

    def _write_table(self, table: pa.Table):
        """Route tables that still carry ``intensities`` through the record path."""
        if "intensities" in table.schema.names:
            self._write_batch(table.to_pylist())
            return
        super()._write_table(table)

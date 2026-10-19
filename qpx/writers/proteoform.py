"""ProteoformWriter — writes proteoform.parquet files.

The proteoform view is flattened like pg: one row per label. Records may carry
the natural ``intensities`` list (``list<{label, intensity}>``); the writer
unrolls it with the same rule as :mod:`qpx.writers.pg`, and an
identification-only member becomes one row with ``label``/``intensity`` null.
"""

from __future__ import annotations

import pyarrow as pa

from qpx.core.data import ProteoformSchema
from qpx.writers.base import BaseWriter
from qpx.writers.pg import _explode_pg_records


class ProteoformWriter(BaseWriter):
    _schema_class = ProteoformSchema

    def _write_batch(self, records: list[dict]):
        super()._write_batch(_explode_pg_records(records))

    def _write_dataframe(self, df):
        if "intensities" in df.columns:
            self._write_batch(df.to_dict("records"))
            return
        super()._write_dataframe(df)

    def _write_table(self, table: pa.Table):
        if "intensities" in table.schema.names:
            exploded = _explode_pg_records(table.to_pylist())
            relaxed = self._relaxed_id_schema()
            batch = pa.RecordBatch.from_pylist(exploded, schema=relaxed)
            super()._write_table(pa.Table.from_batches([batch], schema=relaxed))
            return
        super()._write_table(table)

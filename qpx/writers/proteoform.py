"""ProteoformWriter — writes proteoform.parquet files.

The proteoform view is flattened like pg: one row per label. Records may carry
the natural ``intensities`` list (``list<{label, intensity}>``); the writer
unrolls it with the same rule as :mod:`qpx.writers.pg`, and an
identification-only member becomes one row with ``label``/``intensity`` null.
"""

from qpx.core.data import ProteoformSchema
from qpx.writers.pg import PgWriter


class ProteoformWriter(PgWriter):
    """Write proteoform.parquet files, one row per label as in pg."""

    _schema_class = ProteoformSchema

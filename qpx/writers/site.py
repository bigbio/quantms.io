"""SiteWriter — writes site.parquet files."""

from qpx.core.data import SiteSchema
from qpx.writers.base import BaseWriter


class SiteWriter(BaseWriter):
    """Write schema-validated site.parquet files."""

    _schema_class = SiteSchema

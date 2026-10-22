"""QueryResult — thin wrapper over DuckDB query results for materialization."""

import pandas as pd
import pyarrow as pa


class QueryResult:
    """Wraps a raw DuckDB execution result and provides materialization methods.

    Supports conversion to pandas DataFrame, Arrow Table, Polars DataFrame,
    and row-level iteration via fetchone/fetchall/__iter__.
    """

    def __init__(self, raw_result):
        self._raw = raw_result

    def to_df(self) -> pd.DataFrame:
        """Materialize as a pandas DataFrame."""
        return self._raw.fetchdf()

    def to_arrow(self) -> pa.Table:
        """Materialize as an Arrow Table."""
        return self._raw.to_arrow_table()

    def to_polars(self):
        """Materialize as a Polars DataFrame."""
        import polars as pl

        return pl.from_arrow(self.to_arrow())

    def fetchone(self):
        """Fetch a single row as a tuple."""
        return self._raw.fetchone()

    def fetchall(self):
        """Fetch all rows as a list of tuples."""
        return self._raw.fetchall()

    def __iter__(self):
        """Iterate over rows."""
        while True:
            row = self._raw.fetchone()
            if row is None:
                break
            yield row

"""QueryResult — thin wrapper over DuckDB query results for materialization."""

import pandas as pd
import pyarrow as pa


class QueryResult:
    """Read pandas, Arrow and row output through one DuckDB relation."""

    def __init__(self, connection, sql: str):
        # sql() runs a statement that is not a query at once; DuckDB 1.5 then returns None.
        self._relation = connection.sql(sql)
        self._rows = None

    def to_df(self) -> pd.DataFrame:
        """Materialize as a pandas DataFrame."""
        if self._relation is None:
            return pd.DataFrame()
        return self._relation.df()

    def to_arrow(self) -> pa.Table:
        """Materialize as an Arrow Table."""
        if self._relation is None:
            return pa.table({})
        return self._relation.to_arrow_table()

    def to_polars(self):
        """Materialize as a Polars DataFrame."""
        import polars as pl

        return pl.from_arrow(self.to_arrow())

    def _row_iter(self):
        # Rows are read in one go: a partly read execute() result blocks relation reads on DuckDB 2.0.
        if self._rows is None:
            self._rows = iter([] if self._relation is None else self._relation.fetchall())
        return self._rows

    def fetchone(self):
        """Fetch a single row as a tuple."""
        return next(self._row_iter(), None)

    def fetchall(self):
        """Fetch all rows as a list of tuples."""
        return list(self._row_iter())

    def __iter__(self):
        """Iterate over rows."""
        return self._row_iter()

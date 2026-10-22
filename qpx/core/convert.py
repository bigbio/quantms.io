"""QueryResult — thin wrapper over DuckDB query results for materialization."""

import pandas as pd
import pyarrow as pa


class QueryResult:
    """Read Arrow and pandas output through a relation, and rows through a cursor."""

    def __init__(self, connection, sql: str):
        self._connection = connection
        self._sql = sql
        self._relation = connection.sql(sql)
        self._cursor = None

    def to_df(self) -> pd.DataFrame:
        """Materialize as a pandas DataFrame."""
        return self._relation.df()

    def to_arrow(self) -> pa.Table:
        """Materialize as an Arrow Table."""
        return self._relation.to_arrow_table()

    def to_polars(self):
        """Materialize as a Polars DataFrame."""
        import polars as pl

        return pl.from_arrow(self.to_arrow())

    def _get_cursor(self):
        # DuckDB 2.0 relations have no fetchone().
        if self._cursor is None:
            self._cursor = self._connection.execute(self._sql)
        return self._cursor

    def fetchone(self):
        """Fetch a single row as a tuple."""
        return self._get_cursor().fetchone()

    def fetchall(self):
        """Fetch all rows as a list of tuples."""
        return self._get_cursor().fetchall()

    def __iter__(self):
        """Iterate over rows."""
        cursor = self._get_cursor()
        while True:
            row = cursor.fetchone()
            if row is None:
                break
            yield row

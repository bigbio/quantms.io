"""LazyQuery builder — immutable, chainable SQL query construction."""

from __future__ import annotations

from typing import TYPE_CHECKING

from qpx.core.sql import sql_build, validate_identifier, validate_table

if TYPE_CHECKING:
    from qpx.core.engine import DuckDBEngine


def _escape_sql_string(value: str) -> str:
    """Escape single quotes in a SQL string literal to prevent injection."""
    return value.replace("'", "''")


class LazyQuery:
    """
    Immutable lazy SQL query builder.

    Every filter/select/join returns a NEW LazyQuery — the original is
    unchanged. Nothing hits DuckDB until materialize() is called.

    This is the query engine shared by ALL data structures. Optimizations
    (predicate pushdown hints, query caching, etc.) go here.
    """

    def __init__(self, engine: "DuckDBEngine", table_name: str, sql: str | None = None):
        self._engine = engine
        self._table_name = table_name
        self._sql = sql  # None means "SELECT * FROM table_name"

    @property
    def source(self) -> str:
        """SQL fragment for use in FROM clause."""
        if self._sql:
            return f"({self._sql})"
        return self._table_name

    def select(self, *columns: str) -> LazyQuery:
        cols = ", ".join(columns)
        return LazyQuery(
            self._engine,
            self._table_name,
            sql_build("SELECT $cols FROM $src", cols=cols, src=self.source),
        )

    def filter(self, condition: str) -> LazyQuery:
        return LazyQuery(
            self._engine,
            self._table_name,
            sql_build("SELECT * FROM $src WHERE $cond", src=self.source, cond=condition),
        )

    def limit(self, n: int) -> LazyQuery:
        return LazyQuery(
            self._engine,
            self._table_name,
            sql_build("SELECT * FROM $src LIMIT $n", src=self.source, n=str(int(n))),
        )

    def join(self, other: LazyQuery, on: str) -> LazyQuery:
        stmt = sql_build(
            "SELECT * FROM $src_a a JOIN $src_b b ON a.$col = b.$col",
            src_a=self.source,
            src_b=other.source,
            col=on,
        )
        return LazyQuery(self._engine, self._table_name, stmt)

    def join_list_membership(
        self,
        other: LazyQuery,
        list_column: str,
        value_column: str,
    ) -> LazyQuery:
        """Join rows by expanding a list column and matching each member."""
        stmt = sql_build(
            """
            SELECT a.*, b.*
            FROM $src_a a
            CROSS JOIN UNNEST(a.$list_col) AS _items(item)
            JOIN $src_b b ON _items.item = b.$value_col
            """,
            src_a=self.source,
            src_b=other.source,
            list_col=validate_identifier(list_column),
            value_col=validate_identifier(value_column),
        )
        return LazyQuery(self._engine, self._table_name, stmt)

    def order_by(self, *columns: str, desc: bool = False) -> LazyQuery:
        direction = "DESC" if desc else "ASC"
        cols = ", ".join(c + " " + direction for c in columns)
        return LazyQuery(
            self._engine,
            self._table_name,
            sql_build("SELECT * FROM $src ORDER BY $cols", src=self.source, cols=cols),
        )

    def group_by(self, *columns: str, agg: str = "COUNT(*) AS count") -> LazyQuery:
        cols = ", ".join(columns)
        return LazyQuery(
            self._engine,
            self._table_name,
            sql_build(
                "SELECT $cols, $agg FROM $src GROUP BY $cols",
                cols=cols,
                agg=agg,
                src=self.source,
            ),
        )

    def distinct_values(self, column: str) -> list:
        """Return distinct values for a column."""
        stmt = sql_build(
            "SELECT DISTINCT $col FROM $src ORDER BY $col",
            col=column,
            src=self.source,
        )
        result = self._engine.execute(stmt)
        return [row[0] for row in result.fetchall()]

    def build_sql(self) -> str:
        if self._sql:
            return self._sql
        return sql_build("SELECT * FROM $tbl", tbl=validate_table(self._table_name))

    def execute(self):
        """Execute the query and return raw DuckDB result."""
        return self._engine.execute(self.build_sql())

    def count(self) -> int:
        stmt = sql_build("SELECT COUNT(*) FROM $src", src=self.source)
        row = self._engine.execute(stmt).fetchall()[0]
        return row[0]

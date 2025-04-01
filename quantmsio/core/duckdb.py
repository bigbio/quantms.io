from pathlib import Path
from typing import Union, Generator, Optional

import duckdb
import time
import logging
import os

import pandas as pd

from quantmsio.core.project import create_uuid_filename


class DuckDB:

    def __init__(
        self,
        report_path: Union[Path, str],
        duckdb_max_memory: str = "16GB",
        duckdb_threads: int = 4,
    ):
        self.parquet_db: duckdb.DuckDBPyConnection | None = None
        self._report_path = report_path
        self._duckdb_name = create_uuid_filename("report-duckdb", ".db")
        self._duckdb = self.create_duckdb_from_diann_report(
            duckdb_max_memory, duckdb_threads
        )

    def create_duckdb_from_diann_report(
        self, max_memory: Optional[str], worker_threads: int | None
    ) -> duckdb.DuckDBPyConnection:
        s = time.time()

        database: duckdb.DuckDBPyConnection = duckdb.connect(self._duckdb_name)
        database.execute("SET enable_progress_bar=true")

        if max_memory is not None:
            database.execute("SET max_memory='{}'".format(max_memory))
        if worker_threads is not None:
            database.execute("SET worker_threads='{}'".format(worker_threads))

        msg = database.execute(
            "SELECT * FROM duckdb_settings() where name in ('worker_threads', 'max_memory')"
        ).df()
        logging.info("duckdb uses {} threads.".format(str(msg["value"][0])))
        logging.info("duckdb uses {} of memory.".format(str(msg["value"][1])))

        database.execute(
            "CREATE TABLE report AS SELECT * FROM '{}'".format(self._report_path)
        )
        et = time.time() - s
        logging.info("Time to create duckdb database {} seconds".format(et))
        return database

    def iter_file(
        self, filed: str, file_num: int = 10, columns: list | None = None
    ) -> Generator[tuple[list, pd.DataFrame], None, None]:

        references: list = self.get_unique_references(filed)
        ref_list = [
            references[i : i + file_num] for i in range(0, len(references), file_num)
        ]
        for refs in ref_list:
            batch_df = self.get_report(filed, refs, columns)
            yield refs, batch_df

    def get_unique_references(self, field: str) -> list:

        unique_reference = self._duckdb.sql(f"SELECT DISTINCT {field} FROM report").df()
        return unique_reference[field].tolist()

    def get_report(
        self, filed: str, runs: list, columns: list | None = None
    ) -> pd.DataFrame:

        cols = ", ".join(columns) if columns and isinstance(columns, list) else "*"
        cols = cols.replace("unique", '"unique"')
        database = self.parquet_db.sql(
            """
            select {} from report
            where {} IN {}
            """.format(
                cols, filed, tuple(runs)
            )
        )
        report = database.df()
        return report

    def query_field(
        self, field: str, querys: list, columns: list | None = None
    ) -> pd.DataFrame:

        proteins_key: list = [f"{field} LIKE '%{p}%'" for p in querys]
        query_key: str = " OR ".join(proteins_key)
        cols = ", ".join(columns) if columns and isinstance(columns, list) else "*"
        database = self._duckdb.sql(f"SELECT {cols} FROM report WHERE {query_key}")
        return database.df()

    def destroy_duckdb_database(self) -> None:
        if self._duckdb_name and self._duckdb:
            self._duckdb.close()
            os.remove(self._duckdb_name)
            self._duckdb_name = None
            self._duckdb = None

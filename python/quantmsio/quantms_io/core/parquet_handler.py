import pyarrow as pa
import pyarrow.parquet as pq



class ParquetHandler:
    """
    Base class for handling parquet files
    """

    def __init__(self, parquet_path: str = None):
        self.schema = self._create_schema()
        self.parquet_path = parquet_path
        self.dataset = None

    def _create_schema(self) -> pa.Schema:
        pass

    def write_single_file_parquet(self, table: pa.Table, parquet_output: str = None):
        if self.parquet_path is None and parquet_output is None:
            raise ValueError("parquet_output is required")
        if parquet_output is None:
            parquet_output = self.parquet_path

        pq.write_table(table, parquet_output, compression="snappy")

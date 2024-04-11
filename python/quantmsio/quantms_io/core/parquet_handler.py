import pyarrow as pa
import pyarrow.parquet as pq


class ParquetHandler:
    """
    Base class for handling parquet files
    """

    COMPRESSION = "snappy"  # Compression algorithm for parquet files

    def __init__(self, parquet_path: str = None):
        self.schema = self._create_schema()
        self.parquet_path = parquet_path
        self.dataset = None

    def _create_schema(self) -> pa.Schema:
        """
        Each subclass should implement this method to create the schema for the parquet file.
        """
        pass

    def write_single_file_parquet(self, table: pa.Table, parquet_output: str = None, write_metadata: bool = False):
        """
        Write a single file parquet file. If parquet_output is not specified, the parquet_path will be used. Pyarrow
        is still not writing the metadata of each field in the schema (
        https://stackoverflow.com/questions/76087095/why-pyarrow-write-dataset-not-keeping-field-metadata),
        so we need to allow creating a new parquet file with the metadata.
        Todo: write_metadata should be explored because the current write is not keeping the metadata of each field.
        :param table: pyarrow Table :param
        :param parquet_output: Output parquet file
        :param write_metadata: Write metadata to the parquet file
        """
        if self.parquet_path is None and parquet_output is None:
            raise ValueError("parquet_output is required")
        if parquet_output is None:
            parquet_output = self.parquet_path

        pq.write_table(
            table,
            parquet_output,
            compression=ParquetHandler.COMPRESSION,
            store_schema=True,
        )

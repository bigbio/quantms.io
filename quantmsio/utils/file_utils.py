"""
File utility functions for quantmsio.
This module provides functions for file operations, optimized for performance.
"""

import logging
import os
import pyarrow.parquet as pq
import psutil
import pandas as pd
from pathlib import Path
from typing import List, Tuple, Dict, Iterator

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def extract_protein_list(file: str) -> List[str]:
    """
    Extract a list of proteins from a file.

    Parameters:
    -----------
    file : str
        Path to the file containing protein accessions (one per line)

    Returns:
    --------
    List[str]
        List of protein accessions
    """
    # Use context manager for proper file handling
    with open(file, encoding="utf-8") as f:
        protein_list = [line.strip() for line in f]
    return protein_list


def delete_files_extension(folder: str, extension: str) -> None:
    """
    Delete all files with the given extension in the given folder.

    Parameters:
    -----------
    folder : str
        Folder path
    extension : str
        File extension to match (e.g., '.txt')
    """
    # Use Path for more reliable file operations
    folder_path = Path(folder)
    for file_path in folder_path.glob(f"*{extension}"):
        file_path.unlink()
        logger.info(f"Deleted {file_path}")


def extract_len(file_path: str, header: str) -> Tuple[int, int]:
    """
    Extract the length and position of a section in a tab-delimited file.
    Optimized for performance with proper file handling.

    Parameters:
    -----------
    file_path : str
        Path to the file
    header : str
        Header to search for

    Returns:
    --------
    Tuple[int, int]
        Length of the section and position in the file
    """
    map_tag = {"PSH": "PSM", "PEH": "PEP", "PRH": "PRT"}

    # Check file size first
    if os.stat(file_path).st_size == 0:
        raise ValueError(f"File {file_path} is empty")

    # Use context manager for proper file handling
    with open(file_path, "r") as f:
        pos = 0
        # Find the header
        for line in f:
            if line.split("\t")[0] == header:
                break
            pos = f.tell()

        # Count lines in the section
        file_len = 0
        for line in f:
            if line.split("\t")[0] != map_tag[header]:
                break
            file_len += 1

    return file_len, pos


def load_de_or_ae(path: str) -> Tuple[pd.DataFrame, str]:
    """
    Load differential expression or absolute expression file.
    Optimized to handle comments at the beginning of the file.

    Parameters:
    -----------
    path : str
        Path to the DE or AE file

    Returns:
    --------
    Tuple[pd.DataFrame, str]
        DataFrame containing the data and string containing the comments
    """
    content = []

    # Use context manager for proper file handling
    with open(path, encoding="utf-8") as f:
        # Read comment lines
        for line in f:
            if not line.startswith("#"):
                break
            content.append(line)

        # Reset file position to beginning
        f.seek(0)

        # Skip comment lines when reading with pandas
        df = pd.read_csv(f, sep="\t", comment="#")

    return df, "".join(content)


def read_large_parquet(
    parquet_path: str, batch_size: int = 500000
) -> Iterator[pd.DataFrame]:
    """
    Read a large parquet file in batches to reduce memory usage.

    Parameters:
    -----------
    parquet_path : str
        Path to the parquet file
    batch_size : int, optional
        Number of rows to read in each batch, defaults to 500000

    Yields:
    -------
    pd.DataFrame
        Batch of data from the parquet file
    """
    parquet_file = pq.ParquetFile(parquet_path)

    # Use memory-efficient batch processing
    for batch in parquet_file.iter_batches(batch_size=batch_size):
        yield batch.to_pandas()


def calculate_buffer_size(file_path: str) -> int:
    # Get the total available system memory
    total_memory = psutil.virtual_memory().available

    # Get the size of the file
    file_size = os.path.getsize(file_path)

    # Set the buffer size based on a fraction of available memory or a maximum size
    max_buffer_size = 1 * 1024 * 1024 * 1024  # 1GB
    fraction_of_memory = 0.4  # Adjust as needed

    return min(int(total_memory * fraction_of_memory), max_buffer_size, file_size)


def save_slice_file(
    parquet_table, pqwriters: Dict, output_folder: str, partitions, filename: str
) -> Dict:
    """
    Save a parquet table to a file with partitioning.

    Parameters:
    -----------
    parquet_table : pyarrow.Table
        Table to save
    pqwriters : Dict
        Dictionary of parquet writers
    output_folder : str
        Base output folder
    partitions : tuple
        Partition values
    filename : str
        Output filename

    Returns:
    --------
    Dict
        Updated dictionary of parquet writers
    """
    # Create folder path using Path for better cross-platform compatibility
    folder_path = Path(output_folder)
    for part in partitions:
        folder_path = folder_path / str(part)

    # Create directory if it doesn't exist
    folder_path.mkdir(parents=True, exist_ok=True)

    # Create file path
    save_path = folder_path / filename

    # Create writer if it doesn't exist
    if partitions not in pqwriters:
        pqwriters[partitions] = pq.ParquetWriter(str(save_path), parquet_table.schema)

    # Write table
    pqwriters[partitions].write_table(parquet_table)

    return pqwriters


def save_file(parquet_table, pqwriter, output_folder: str, filename: str):
    """
    Save a parquet table to a file.

    Parameters:
    -----------
    parquet_table : pyarrow.Table
        Table to save
    pqwriter : ParquetWriter or None
        Existing parquet writer or None
    output_folder : str
        Output folder
    filename : str
        Output filename

    Returns:
    --------
    ParquetWriter
        Parquet writer object
    """
    # Create directory if it doesn't exist
    folder_path = Path(output_folder)
    folder_path.mkdir(parents=True, exist_ok=True)

    # Create file path
    save_path = folder_path / filename

    # Create writer if it doesn't exist
    if not pqwriter:
        pqwriter = pq.ParquetWriter(str(save_path), parquet_table.schema)

    # Write table
    pqwriter.write_table(parquet_table)

    return pqwriter


def close_file(pqwriters: Dict = None, pqwriter=None) -> None:
    """
    Close parquet writers.

    Parameters:
    -----------
    pqwriters : Dict, optional
        Dictionary of parquet writers
    pqwriter : ParquetWriter, optional
        Single parquet writer
    """
    if pqwriter:
        pqwriter.close()
    elif pqwriters:
        for writer in pqwriters.values():
            writer.close()


def find_ae_files(directory: str) -> List[Path]:
    """
    Find absolute expression files in a directory.

    Parameters:
    -----------
    directory : str
        Directory to search

    Returns:
    --------
    List[Path]
        List of absolute expression file paths
    """
    path = Path(directory)
    ae_files = list(path.rglob("*.absolute.tsv"))
    return ae_files

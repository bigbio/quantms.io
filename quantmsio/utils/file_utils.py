import logging
import os
import pyarrow.parquet as pq
import psutil
import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def extract_protein_list(file):
    f = open(file, encoding="utf-8")
    protein_list = f.read().splitlines()
    return protein_list


def delete_files_extension(folder: str, extension: str) -> None:
    """
    Delete all files with the given extension in the given folder
    :param folder: Folder path
    :param extension: Extension
    """
    for file in os.listdir(folder):
        if file.endswith(extension):
            os.remove(f"{folder}/{file}")
            logger.info(f"Deleted {folder}/{file}")


def extract_len(fle, header):
    map_tag = {"PSH": "PSM", "PEH": "PEP", "PRH": "PRT"}
    if os.stat(fle).st_size == 0:
        raise ValueError("File is empty")
    f = open(fle)
    pos = 0
    line = f.readline()
    while line.split("\t")[0] != header:
        pos = f.tell()
        line = f.readline()
    line = f.readline()
    fle_len = 0
    while line.split("\t")[0] == map_tag[header]:
        fle_len += 1
        line = f.readline()
    f.close()
    return fle_len, pos


def load_de_or_ae(path):
    f = open(path, encoding="utf-8")
    line = f.readline()
    pos = 0
    content = ""
    while line.startswith("#"):
        pos = f.tell()
        content += line
        line = f.readline()
    f.seek(pos - 1)
    return pd.read_csv(f, sep="\t"), content

def read_large_parquet(parquet_path: str, batch_size: int = 500000):
    """_summary_
    :param parquet_path: _description_
    :param batch_size: _description_, defaults to 100000
    :yield: _description_
    """
    parquet_file = pq.ParquetFile(parquet_path)
    for batch in parquet_file.iter_batches(batch_size=batch_size):
        batch_df = batch.to_pandas()
        yield batch_df

def calculate_buffer_size(file_path: str) -> int:
    # Get the total available system memory
    total_memory = psutil.virtual_memory().available

    # Get the size of the file
    file_size = os.path.getsize(file_path)

    # Set the buffer size based on a fraction of available memory or a maximum size
    max_buffer_size = 1 * 1024 * 1024 * 1024  # 1GB
    fraction_of_memory = 0.4  # Adjust as needed

    return min(int(total_memory * fraction_of_memory), max_buffer_size, file_size)

import os

import logging

import psutil

logging.basicConfig(level = logging.INFO)
logger = logging.getLogger(__name__)

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

def calculate_buffer_size(file_path: str) -> int:
    # Get the total available system memory
    total_memory = psutil.virtual_memory().available

    # Get the size of the file
    file_size = os.path.getsize(file_path)

    # Set the buffer size based on a fraction of available memory or a maximum size
    max_buffer_size = 1 * 1024 * 1024 * 1024  # 1GB
    fraction_of_memory = 0.4  # Adjust as needed

    return min(int(total_memory * fraction_of_memory), max_buffer_size, file_size)
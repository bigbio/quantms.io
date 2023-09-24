import os

import logging
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
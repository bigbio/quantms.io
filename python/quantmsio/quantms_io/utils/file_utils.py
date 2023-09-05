import os


def delete_files_extension(folder: str, extension: str) -> None:
    """
    Delete all files with the given extension in the given folder
    :param folder: Folder path
    :param extension: Extension
    """
    for file in os.listdir(folder):
        if file.endswith(extension):
            os.remove(f"{folder}/{file}")
            print(f"Deleted {folder}/{file}")

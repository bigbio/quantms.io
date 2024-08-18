from pathlib import Path


TEST_DATA_ROOT = Path(__file__).parent / "examples"


def datafile(path: str):
    path = str(path).removeprefix("/").removeprefix("examples/")
    return str(TEST_DATA_ROOT / path)
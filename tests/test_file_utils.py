import pytest
import tempfile
import os
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path
from unittest.mock import patch, MagicMock
from quantmsio.utils.file_utils import (
    extract_protein_list,
    delete_files_extension,
    extract_len,
    load_de_or_ae,
    read_large_parquet,
    calculate_buffer_size,
    save_slice_file,
    save_file,
    close_file,
    find_ae_files,
)


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as temp_dir:
        yield Path(temp_dir)


def test_extract_protein_list(temp_dir):
    """Test extracting protein list from a file."""
    # Create a test protein list file
    protein_file = temp_dir / "proteins.txt"
    with open(protein_file, "w") as f:
        f.write("Protein1\nProtein2\nProtein3\n")

    # Test the function
    proteins = extract_protein_list(str(protein_file))
    assert len(proteins) == 3
    assert proteins == ["Protein1", "Protein2", "Protein3"]


def test_delete_files_extension(temp_dir):
    """Test deleting files with a specific extension."""
    # Create test files
    (temp_dir / "file1.txt").touch()
    (temp_dir / "file2.txt").touch()
    (temp_dir / "file3.csv").touch()

    # Test the function
    delete_files_extension(str(temp_dir), ".txt")

    # Check that only .txt files were deleted
    files = list(temp_dir.glob("*"))
    assert len(files) == 1
    assert files[0].name == "file3.csv"


def test_extract_len(temp_dir):
    """Test extracting length and position from a tab-delimited file."""
    # Create a test file
    test_file = temp_dir / "test.txt"
    with open(test_file, "w") as f:
        f.write("Header\tValue\n")
        f.write("PSH\tHeader\n")
        f.write("PSM\tValue1\n")
        f.write("PSM\tValue2\n")
        f.write("PSM\tValue3\n")
        f.write("Other\tValue\n")

    # Test the function
    length, pos = extract_len(str(test_file), "PSH")
    assert length == 3  # 3 PSM lines


def test_load_de_or_ae(temp_dir):
    """Test loading differential expression or absolute expression file."""
    # Create a test file
    test_file = temp_dir / "test.tsv"
    with open(test_file, "w") as f:
        f.write("# Comment 1\n")
        f.write("# Comment 2\n")
        f.write("Column1\tColumn2\n")
        f.write("Value1\tValue2\n")
        f.write("Value3\tValue4\n")

    # Test the function
    df, comments = load_de_or_ae(str(test_file))
    assert len(df) == 2  # 2 data rows
    assert df.shape[1] == 2  # 2 columns
    assert "# Comment 1" in comments
    assert "# Comment 2" in comments


@patch("pyarrow.parquet.ParquetFile")
def test_read_large_parquet(mock_parquet_file):
    """Test reading a large parquet file in batches."""
    # Mock the ParquetFile and iter_batches
    mock_batch1 = MagicMock()
    mock_batch1.to_pandas.return_value = pd.DataFrame({"col1": [1, 2]})
    mock_batch2 = MagicMock()
    mock_batch2.to_pandas.return_value = pd.DataFrame({"col1": [3, 4]})

    mock_parquet_file.return_value.iter_batches.return_value = [mock_batch1, mock_batch2]

    # Test the function
    batches = list(read_large_parquet("dummy.parquet", 1000))
    assert len(batches) == 2
    assert batches[0]["col1"].tolist() == [1, 2]
    assert batches[1]["col1"].tolist() == [3, 4]


@patch("psutil.virtual_memory")
@patch("os.path.getsize")
def test_calculate_buffer_size(mock_getsize, mock_virtual_memory):
    """Test calculating buffer size based on system memory."""
    # Mock the system memory and file size
    mock_virtual_memory.return_value.available = 8 * 1024 * 1024 * 1024  # 8GB
    mock_getsize.return_value = 500 * 1024 * 1024  # 500MB

    # Test the function
    buffer_size = calculate_buffer_size("dummy.file")

    # Should be 40% of available memory, but capped at 1GB or file size
    assert buffer_size == 500 * 1024 * 1024  # Should be file size (500MB)


def test_save_slice_file_and_close_file(temp_dir):
    """Test saving a parquet table with partitioning and closing writers."""
    # Create test data
    data = {"col1": [1, 2, 3], "col2": ["a", "b", "c"]}
    df = pd.DataFrame(data)
    table = pa.Table.from_pandas(df)

    # Mock ParquetWriter
    mock_writer = MagicMock()
    pqwriters = {}

    with patch("pyarrow.parquet.ParquetWriter", return_value=mock_writer):
        # Test save_slice_file
        pqwriters = save_slice_file(table, pqwriters, str(temp_dir), ("part1", "part2"), "test.parquet")

        # Check that the writer was created and write_table was called
        assert ("part1", "part2") in pqwriters
        mock_writer.write_table.assert_called_once_with(table)

        # Test close_file
        close_file(pqwriters=pqwriters)
        mock_writer.close.assert_called_once()


def test_save_file(temp_dir):
    """Test saving a parquet table to a file."""
    # Create test data
    data = {"col1": [1, 2, 3], "col2": ["a", "b", "c"]}
    df = pd.DataFrame(data)
    table = pa.Table.from_pandas(df)

    # Mock ParquetWriter
    mock_writer = MagicMock()

    with patch("pyarrow.parquet.ParquetWriter", return_value=mock_writer):
        # Test with no existing writer
        writer = save_file(table, None, str(temp_dir), "test.parquet")

        # Check that the writer was created and write_table was called
        assert writer == mock_writer
        mock_writer.write_table.assert_called_once_with(table)

        # Test with existing writer
        mock_writer.reset_mock()
        writer = save_file(table, mock_writer, str(temp_dir), "test.parquet")

        # Check that write_table was called but no new writer was created
        assert writer == mock_writer
        mock_writer.write_table.assert_called_once_with(table)


def test_find_ae_files(temp_dir):
    """Test finding absolute expression files in a directory."""
    # Create test files
    (temp_dir / "file1.absolute.tsv").touch()
    (temp_dir / "file2.absolute.tsv").touch()
    (temp_dir / "file3.tsv").touch()
    (temp_dir / "subdir").mkdir()
    (temp_dir / "subdir" / "file4.absolute.tsv").touch()

    # Test the function
    ae_files = find_ae_files(str(temp_dir))

    # Should find 3 .absolute.tsv files (including in subdirectory)
    assert len(ae_files) == 3
    file_names = [f.name for f in ae_files]
    assert "file1.absolute.tsv" in file_names
    assert "file2.absolute.tsv" in file_names
    assert "file4.absolute.tsv" in file_names

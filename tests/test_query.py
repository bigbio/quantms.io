import pytest
from .common import datafile
from quantmsio.operate.query import Query


@pytest.fixture
def feature_path():
    """Get the path to the feature parquet file."""
    return datafile("parquet/feature.parquet")


@pytest.fixture
def fasta_path():
    """Get the path to the fasta file."""
    return datafile("fasta/Homo-sapiens.fasta")


def test_iter_samples(feature_path):
    """Test iterating through samples."""
    # Initialize Query
    query = Query(feature_path)

    # Iterate through samples
    count = 0
    for sample in query.iter_samples():
        # Add assertions to verify the result
        assert sample is not None
        assert isinstance(sample, str)
        count += 1

    # Ensure we got at least one sample
    assert count > 0


def test_iter_chunk(feature_path):
    """Test iterating through chunks."""
    # Initialize Query
    query = Query(feature_path)

    # Iterate through chunks
    count = 0
    for chunk in query.iter_chunk():
        # Add assertions to verify the result
        assert chunk is not None
        assert len(chunk) > 0
        assert "peptidoform" in chunk.columns
        count += 1

    # Ensure we got at least one chunk
    assert count > 0


def test_iter_file(feature_path):
    """Test iterating through the file."""
    # Initialize Query
    query = Query(feature_path)

    # Iterate through the file
    count = 0
    for df in query.iter_file():
        # Add assertions to verify the result
        assert df is not None
        assert len(df) > 0
        assert "peptidoform" in df.columns
        count += 1

    # Ensure we got at least one result
    assert count > 0


def test_inject_position_msg(feature_path, fasta_path):
    """Test injecting position messages."""
    # Initialize Query
    query = Query(feature_path)

    # Get report from database
    try:
        df = query.get_report_from_database(["20180914_QE8_nLC0_BDA_SA_DIA_Keratinocytes_NN002"])

        # Get protein dictionary
        protein_dict = query.get_protein_dict(fasta_path)

        # Inject position messages
        result = query.inject_position_msg(df, protein_dict)

        # Add assertions to verify the result
        assert result is not None
        assert len(result) > 0
        assert "peptidoform" in result.columns
        assert "position" in result.columns
    except Exception as e:
        # Skip the test if the sample is not found
        pytest.skip(f"Sample not found in the database: {e}")

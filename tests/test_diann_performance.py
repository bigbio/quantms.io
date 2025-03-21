import pytest
import pandas as pd
import numpy as np
import os
import tempfile
import time
from pathlib import Path
from unittest.mock import patch, MagicMock
from quantmsio.core.diann import DiaNNConvert


@pytest.fixture
def mock_setup():
    """Set up test data and mocks."""
    # Create a mock DuckDB connection
    with patch("quantmsio.core.diann.DuckDB.__init__") as mock_duckdb, \
         patch.object(DiaNNConvert, "get_report_from_database") as mock_query, \
         patch.object(DiaNNConvert, "get_masses_and_modifications_map") as mock_masses, \
         patch.object(DiaNNConvert, "add_additional_msg") as mock_add_msg:
        
        mock_duckdb.return_value = None
        
        # Create sample data for testing
        sample_report = pd.DataFrame(
            {
                "run": np.repeat(["run1", "run2"], 500),
                "peptidoform": [f"PEPTIDE{i}" for i in range(1000)],
                "precursor_charge": np.random.choice([2, 3, 4], 1000),
                "qvalue": np.random.random(1000) * 0.1,
                "pg_accessions": [f"PROTEIN{i % 100}" for i in range(1000)],
                "rt": np.random.random(1000) * 100,
            }
        )
        
        # Set up mock return values
        mock_query.return_value = sample_report
        mock_masses.return_value = (
            {f"PEPTIDE{i}": 1000 + i for i in range(1000)},  # masses_map
            {f"PEPTIDE{i}": f"PEPTIDE{i}_mod" for i in range(1000)},  # modifications_map
        )
        
        # Create a temporary directory for MS info files
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create mock MS info files
            for run in ["run1", "run2"]:
                ms_info = pd.DataFrame(
                    {
                        "rt": np.random.random(100) * 6000,  # in seconds
                        "scan": [f"scan_{i}" for i in range(100)],
                        "observed_mz": np.random.random(100) * 1000,
                    }
                )
                ms_info_path = Path(temp_dir) / f"{run}_ms_info.parquet"
                ms_info.to_parquet(ms_info_path)
            
            # Initialize the DiaNNConvert object
            diann = DiaNNConvert("dummy_report.tsv")
            
            yield {
                "diann": diann,
                "temp_dir": temp_dir,
                "sample_report": sample_report,
                "mock_add_msg": mock_add_msg
            }


def test_main_report_df_performance(mock_setup):
    """Test the performance of main_report_df method."""
    diann = mock_setup["diann"]
    temp_dir = mock_setup["temp_dir"]
    sample_report = mock_setup["sample_report"]
    
    # Patch os.listdir to return our mock files
    with patch("os.listdir", return_value=[f"{run}_ms_info.parquet" for run in ["run1", "run2"]]):
        # Measure execution time
        start_time = time.time()
        reports = list(diann.main_report_df(0.05, temp_dir, 2))
        execution_time = time.time() - start_time
        
        # Print execution time for reference
        print(f"main_report_df execution time: {execution_time:.6f} seconds")
        
        # Verify the results
        assert len(reports) == 1  # Should have one batch of reports
        report = reports[0]
        
        # Check that the report has the expected columns
        expected_columns = list(sample_report.columns) + ["scan", "observed_mz", "calculated_mz"]
        for col in expected_columns:
            assert col in report.columns
        
        # Check that peptidoform values were mapped correctly
        assert all(pep.endswith("_mod") for pep in report["peptidoform"])
        
        # Check that calculated_mz was computed
        assert not report["calculated_mz"].isna().any()


def test_generate_feature_performance(mock_setup):
    """Test the performance of generate_feature method."""
    diann = mock_setup["diann"]
    temp_dir = mock_setup["temp_dir"]
    mock_add_msg = mock_setup["mock_add_msg"]
    
    # Patch os.listdir to return our mock files
    with patch("os.listdir", return_value=[f"{run}_ms_info.parquet" for run in ["run1", "run2"]]):
        # Measure execution time
        start_time = time.time()
        features = list(diann.generate_feature(0.05, temp_dir, 2))
        execution_time = time.time() - start_time
        
        # Print execution time for reference
        print(f"generate_feature execution time: {execution_time:.6f} seconds")
        
        # Verify the results
        assert len(features) == 1  # Should have one batch of features
        
        # Check that add_additional_msg was called
        mock_add_msg.assert_called_once()

import unittest
import pandas as pd
import numpy as np
import os
import tempfile
import time
from pathlib import Path
from unittest.mock import patch, MagicMock
from quantmsio.core.diann import DiaNNConvert


class TestDiaNNPerformance(unittest.TestCase):
    """Test performance optimizations in the DiaNNConvert class."""

    def setUp(self):
        """Set up test data and mocks."""
        # Create a mock DuckDB connection
        self.mock_duckdb_patcher = patch('quantmsio.core.diann.DuckDB.__init__')
        self.mock_duckdb = self.mock_duckdb_patcher.start()
        self.mock_duckdb.return_value = None
        
        # Create a mock for the query method
        self.mock_query_patcher = patch.object(DiaNNConvert, 'get_report_from_database')
        self.mock_query = self.mock_query_patcher.start()
        
        # Create a mock for the get_masses_and_modifications_map method
        self.mock_masses_patcher = patch.object(DiaNNConvert, 'get_masses_and_modifications_map')
        self.mock_masses = self.mock_masses_patcher.start()
        
        # Create sample data for testing
        self.sample_report = pd.DataFrame({
            'run': np.repeat(['run1', 'run2'], 500),
            'peptidoform': [f'PEPTIDE{i}' for i in range(1000)],
            'precursor_charge': np.random.choice([2, 3, 4], 1000),
            'qvalue': np.random.random(1000) * 0.1,
            'pg_accessions': [f'PROTEIN{i % 100}' for i in range(1000)],
            'rt': np.random.random(1000) * 100,
        })
        
        # Set up mock return values
        self.mock_query.return_value = self.sample_report
        self.mock_masses.return_value = (
            {f'PEPTIDE{i}': 1000 + i for i in range(1000)},  # masses_map
            {f'PEPTIDE{i}': f'PEPTIDE{i}_mod' for i in range(1000)}  # modifications_map
        )
        
        # Create a temporary directory for MS info files
        self.temp_dir = tempfile.TemporaryDirectory()
        
        # Create mock MS info files
        for run in ['run1', 'run2']:
            ms_info = pd.DataFrame({
                'rt': np.random.random(100) * 6000,  # in seconds
                'scan': [f'scan_{i}' for i in range(100)],
                'observed_mz': np.random.random(100) * 1000,
            })
            ms_info_path = Path(self.temp_dir.name) / f"{run}_ms_info.parquet"
            ms_info.to_parquet(ms_info_path)
        
        # Initialize the DiaNNConvert object
        self.diann = DiaNNConvert('dummy_report.tsv')
        
        # Mock the add_additional_msg method to avoid needing SDRF data
        self.mock_add_msg_patcher = patch.object(DiaNNConvert, 'add_additional_msg')
        self.mock_add_msg = self.mock_add_msg_patcher.start()

    def tearDown(self):
        """Clean up after tests."""
        self.mock_duckdb_patcher.stop()
        self.mock_query_patcher.stop()
        self.mock_masses_patcher.stop()
        self.mock_add_msg_patcher.stop()
        self.temp_dir.cleanup()

    def test_main_report_df_performance(self):
        """Test the performance of main_report_df method."""
        # Patch os.listdir to return our mock files
        with patch('os.listdir', return_value=[f"{run}_ms_info.parquet" for run in ['run1', 'run2']]):
            # Measure execution time
            start_time = time.time()
            reports = list(self.diann.main_report_df(0.05, self.temp_dir.name, 2))
            execution_time = time.time() - start_time
            
            # Print execution time for reference
            print(f"main_report_df execution time: {execution_time:.6f} seconds")
            
            # Verify the results
            self.assertEqual(len(reports), 1)  # Should have one batch of reports
            report = reports[0]
            
            # Check that the report has the expected columns
            expected_columns = list(self.sample_report.columns) + ['scan', 'observed_mz', 'calculated_mz']
            for col in expected_columns:
                self.assertIn(col, report.columns)
            
            # Check that peptidoform values were mapped correctly
            self.assertTrue(all(pep.endswith('_mod') for pep in report['peptidoform']))
            
            # Check that calculated_mz was computed
            self.assertFalse(report['calculated_mz'].isna().any())

    def test_generate_feature_performance(self):
        """Test the performance of generate_feature method."""
        # Patch os.listdir to return our mock files
        with patch('os.listdir', return_value=[f"{run}_ms_info.parquet" for run in ['run1', 'run2']]):
            # Measure execution time
            start_time = time.time()
            features = list(self.diann.generate_feature(0.05, self.temp_dir.name, 2))
            execution_time = time.time() - start_time
            
            # Print execution time for reference
            print(f"generate_feature execution time: {execution_time:.6f} seconds")
            
            # Verify the results
            self.assertEqual(len(features), 1)  # Should have one batch of features
            
            # Check that add_additional_msg was called
            self.mock_add_msg.assert_called_once()


if __name__ == '__main__':
    unittest.main()
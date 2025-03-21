import unittest
import pandas as pd
import numpy as np
import time
from quantmsio.core.feature import Feature


class TestFeaturePerformance(unittest.TestCase):
    """Test performance optimizations in the Feature class."""

    def setUp(self):
        """Set up test data."""
        # Create a sample DataFrame with typical feature data
        self.test_df = pd.DataFrame(
            {
                "pg_global_qvalue": np.random.random(1000),
                "unique": np.random.choice([0, 1], 1000),
                "precursor_charge": np.random.choice([2, 3, 4], 1000),
                "calculated_mz": np.random.random(1000) * 1000,
                "observed_mz": np.random.random(1000) * 1000,
                "posterior_error_probability": np.random.random(1000),
                "is_decoy": np.random.choice([0, 1], 1000),
                "scan": [f"scan_{i}" for i in range(1000)],
                "scan_reference_file_name": [f"file_{i % 10}" for i in range(1000)],
                "rt": np.random.random(1000) * 100,
            }
        )

        # Add some NaN values to test robustness
        for col in self.test_df.columns:
            mask = np.random.choice([True, False], 1000, p=[0.05, 0.95])
            self.test_df.loc[mask, col] = np.nan

    def test_convert_to_parquet_format_performance(self):
        """Test the performance of convert_to_parquet_format method."""
        # Make a copy of the test DataFrame
        df_copy = self.test_df.copy()

        # Measure execution time
        start_time = time.time()
        Feature.convert_to_parquet_format(df_copy)
        execution_time = time.time() - start_time

        # Print execution time for reference
        print(f"convert_to_parquet_format execution time: {execution_time:.6f} seconds")

        # Verify the conversion was done correctly
        self.assertEqual(df_copy["unique"].dtype, pd.Int32Dtype())
        self.assertEqual(df_copy["precursor_charge"].dtype, pd.Int32Dtype())
        self.assertEqual(df_copy["is_decoy"].dtype, pd.Int32Dtype())
        self.assertTrue(pd.api.types.is_float_dtype(df_copy["pg_global_qvalue"]))
        self.assertTrue(pd.api.types.is_float_dtype(df_copy["calculated_mz"]))
        self.assertTrue(pd.api.types.is_float_dtype(df_copy["observed_mz"]))
        self.assertTrue(pd.api.types.is_float_dtype(df_copy["posterior_error_probability"]))
        self.assertTrue(pd.api.types.is_float_dtype(df_copy["rt"]))
        self.assertTrue(pd.api.types.is_string_dtype(df_copy["scan"]))
        self.assertTrue(pd.api.types.is_string_dtype(df_copy["scan_reference_file_name"]))

    def test_slice_method(self):
        """Test the slice method."""
        # Add a partition column
        self.test_df["partition_col"] = np.random.choice(["A", "B", "C"], 1000)

        # Test slicing by partition_col
        slices = list(Feature.slice(self.test_df, ["partition_col"]))

        # Verify we have the correct number of slices
        unique_values = self.test_df["partition_col"].unique()
        self.assertEqual(len(slices), len(unique_values))

        # Verify each slice has the correct data
        for key, df_slice in slices:
            self.assertEqual(df_slice["partition_col"].unique()[0], key)
            self.assertEqual(len(df_slice), len(self.test_df[self.test_df["partition_col"] == key]))

    def test_slice_method_multiple_partitions(self):
        """Test the slice method with multiple partition columns."""
        # Add partition columns
        self.test_df["partition_col1"] = np.random.choice(["A", "B"], 1000)
        self.test_df["partition_col2"] = np.random.choice([1, 2], 1000)

        # Test slicing by both partition columns
        slices = list(Feature.slice(self.test_df, ["partition_col1", "partition_col2"]))

        # Verify each slice has the correct data
        for key, df_slice in slices:
            # key should be a tuple of (partition_col1, partition_col2)
            self.assertEqual(len(key), 2)
            self.assertEqual(df_slice["partition_col1"].unique()[0], key[0])
            self.assertEqual(df_slice["partition_col2"].unique()[0], key[1])


if __name__ == "__main__":
    unittest.main()

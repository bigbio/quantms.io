import pytest
import pandas as pd
import numpy as np
import time
from quantmsio.core.feature import Feature


@pytest.fixture
def test_df():
    """Create a sample DataFrame with typical feature data."""
    # Create a sample DataFrame with typical feature data
    df = pd.DataFrame(
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
    for col in df.columns:
        mask = np.random.choice([True, False], 1000, p=[0.05, 0.95])
        df.loc[mask, col] = np.nan
    
    return df


def test_convert_to_parquet_format_performance(test_df):
    """Test the performance of convert_to_parquet_format method."""
    # Make a copy of the test DataFrame
    df_copy = test_df.copy()

    # Measure execution time
    start_time = time.time()
    Feature.convert_to_parquet_format(df_copy)
    execution_time = time.time() - start_time

    # Print execution time for reference
    print(f"convert_to_parquet_format execution time: {execution_time:.6f} seconds")

    # Verify the conversion was done correctly
    assert df_copy["unique"].dtype == pd.Int32Dtype()
    assert df_copy["precursor_charge"].dtype == pd.Int32Dtype()
    assert df_copy["is_decoy"].dtype == pd.Int32Dtype()
    assert pd.api.types.is_float_dtype(df_copy["pg_global_qvalue"])
    assert pd.api.types.is_float_dtype(df_copy["calculated_mz"])
    assert pd.api.types.is_float_dtype(df_copy["observed_mz"])
    assert pd.api.types.is_float_dtype(df_copy["posterior_error_probability"])
    assert pd.api.types.is_float_dtype(df_copy["rt"])
    assert pd.api.types.is_string_dtype(df_copy["scan"])
    assert pd.api.types.is_string_dtype(df_copy["scan_reference_file_name"])


def test_slice_method(test_df):
    """Test the slice method."""
    # Add a partition column
    test_df["partition_col"] = np.random.choice(["A", "B", "C"], 1000)

    # Test slicing by partition_col
    slices = list(Feature.slice(test_df, ["partition_col"]))

    # Verify we have the correct number of slices
    unique_values = test_df["partition_col"].unique()
    assert len(slices) == len(unique_values)

    # Verify each slice has the correct data
    for key, df_slice in slices:
        assert df_slice["partition_col"].unique()[0] == key
        assert len(df_slice) == len(test_df[test_df["partition_col"] == key])


def test_slice_method_multiple_partitions(test_df):
    """Test the slice method with multiple partition columns."""
    # Add partition columns
    test_df["partition_col1"] = np.random.choice(["A", "B"], 1000)
    test_df["partition_col2"] = np.random.choice([1, 2], 1000)

    # Test slicing by both partition columns
    slices = list(Feature.slice(test_df, ["partition_col1", "partition_col2"]))

    # Verify each slice has the correct data
    for key, df_slice in slices:
        # key should be a tuple of (partition_col1, partition_col2)
        assert len(key) == 2
        assert df_slice["partition_col1"].unique()[0] == key[0]
        assert df_slice["partition_col2"].unique()[0] == key[1]

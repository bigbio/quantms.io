from unittest import TestCase
from quantms_io.core.statistics import ParquetStatistics


class TestStatistics(TestCase):

    def test_get_number_of_peptides(self):
        # Mock the duckdb.sql method to return a DataFrame with a single column 'sequence'
        parquet_path = "/Users/yperez/work/quantms-data/PXD037340.2/PXD037340.2-93425767-e297-4e8e-b8dd-30793e3aae19.feature.parquet"
        # Create an instance of ParquetStatistics
        parquet_stats = ParquetStatistics(parquet_path)

        # Call the method get_number_of_peptides
        peptide_count = parquet_stats.get_number_of_peptides()
        peptidoform_count = parquet_stats.get_number_of_peptidoforms()
        sample_accession_count = parquet_stats.get_number_of_samples()
        proteins_count = parquet_stats.get_number_proteins()


        # Assert that the result is equal to the number of unique peptides in the DataFrame
        assert peptide_count == 37783
        assert peptidoform_count == 37783
        assert sample_accession_count == 31
        assert proteins_count == 5578

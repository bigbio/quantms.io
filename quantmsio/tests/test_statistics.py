from unittest import TestCase

from core.statistics import IbaqStatistics
from core.statistics import ParquetStatistics


class TestStatistics(TestCase):

    def test_feature_get_number_of_peptides(self):
        # Mock the duckdb.sql method to return a DataFrame with a single column 'sequence'
        parquet_path = __package__ + "/examples/JSON/PXD040438.feature.parquet"
        # Create an instance of ParquetStatistics
        parquet_stats = ParquetStatistics(parquet_path)

        # Call the method get_number_of_peptides
        peptide_count = parquet_stats.get_number_of_peptides()
        peptidoform_count = parquet_stats.get_number_of_peptidoforms()
        sample_accession_count = parquet_stats.get_number_of_samples()
        proteins_count = parquet_stats.get_number_of_proteins()
        msruns_count = parquet_stats.get_number_msruns()

        # Assert that the result is equal to the number of unique peptides in the DataFrame
        assert peptide_count == 2685
        assert peptidoform_count == 2729
        assert sample_accession_count == 2
        assert proteins_count == 427
        assert msruns_count == 2

    def test_ibaq_get_number_of_proteins(self):
        # Mock the duckdb.sql method to return a DataFrame with a single column 'sequence'
        ibaq_path = __package__ + "/examples/JSON/PXD016999.absolute.tsv"
        # Create an instance of ParquetStatistics
        ibaq_stats = IbaqStatistics(ibaq_path)

        # Call the method get_number_of_peptides
        protein_count = ibaq_stats.get_number_of_proteins()
        sample_accession_count = ibaq_stats.get_number_of_samples()

        # Assert that the result is equal to the number of unique peptides in the DataFrame
        assert protein_count == 1062
        assert sample_accession_count == 1

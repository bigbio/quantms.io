from unittest import TestCase
from quantms_io.core.statistics import ParquetStatistics, IbaqStatistics


class TestStatistics(TestCase):

    def test_feature_get_number_of_peptides(self):
        # Mock the duckdb.sql method to return a DataFrame with a single column 'sequence'
        parquet_path = "/Users/yperez/work/quantms-data/PXD037340.2/PXD037340.2-93425767-e297-4e8e-b8dd-30793e3aae19.feature.parquet"
        # Create an instance of ParquetStatistics
        parquet_stats = ParquetStatistics(parquet_path)

        # Call the method get_number_of_peptides
        peptide_count = parquet_stats.get_number_of_peptides()
        peptidoform_count = parquet_stats.get_number_of_peptidoforms()
        sample_accession_count = parquet_stats.get_number_of_samples()
        proteins_count = parquet_stats.get_number_proteins()
        msruns_count = parquet_stats.get_number_msruns()


        # Assert that the result is equal to the number of unique peptides in the DataFrame
        assert peptide_count == 37783
        assert peptidoform_count == 37783
        assert sample_accession_count == 31
        assert proteins_count == 5578
        assert msruns_count == 110

    def test_ibaq_get_number_of_proteins(self):
        # Mock the duckdb.sql method to return a DataFrame with a single column 'sequence'
        ibaq_path = "/Users/yperez/work/PXD000561-protein-ibaq.tsv"
        # Create an instance of ParquetStatistics
        ibaq_stats = IbaqStatistics(ibaq_path)

        # Call the method get_number_of_peptides
        protein_count = ibaq_stats.get_number_of_proteins()
        sample_accession_count = ibaq_stats.get_number_of_samples()

        # Assert that the result is equal to the number of unique peptides in the DataFrame
        assert protein_count == 12025
        assert sample_accession_count == 85

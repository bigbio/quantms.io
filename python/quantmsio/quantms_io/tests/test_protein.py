from unittest import TestCase, mock

import pandas as pd

from quantms_io.core.protein import ProteinHandler
import pyarrow as pa
import random

class TestProteinHandler(TestCase):
    def test_describe_schema(self):
        protein_handler = ProteinHandler()
        print(protein_handler.describe_schema())

    @mock.patch('pyarrow.parquet.read_table')
    def test_read_proteins(self, mock_read_table):
        # Mock the behavior of pyarrow.parquet.read_table
        mock_table = mock.Mock(spec=pa.Table)
        mock_read_table.return_value = mock_table

        # Create a mock schema
        protein_mock = ProteinHandler()
        schema = protein_mock._create_schema()

        # Create a ProteinHandler instance
        parquet_path = 'mock_protein_data.parquet'
        protein_handler = ProteinHandler(parquet_path)
        protein_handler.schema = schema

        # Test the read_protein_dataset method
        df = protein_handler.read_protein_dataset()

        # Verify that the mock read_table function was called
        mock_read_table.assert_called_once_with(parquet_path)

        # Verify that the returned DataFrame matches the mock_table data
        expected_data = {
            "protein_accessions": [["P123"], ["P456"]],
            "sample_accession": ["S1", "S2"],
            "abundance": [0.5, 0.8],
            # ... other fields ...
        }
        expected_df = pd.DataFrame(expected_data)
        pd.testing.assert_frame_equal(df, expected_df)


    def test_write_an_example_protein(self):
        protein_manager = ProteinHandler()
        protein_manager.parquet_path = 'example_protein.parquet'
        protein_manager.create_protein({
            "protein_accessions": ["P12345"],
            "sample_accession": "S1",
            "abundance": 0.75,
            "global_qvalue": 0.05,
            "is_decoy": 0,
            "best_id_score": "Score123: 0.75",
            "gene_accessions": ["G123", "G456"],
            "gene_names": ["GeneA", "GeneB"],
            "number_of_peptides": 15,
            "number_of_psms": 25,
            "number_of_unique_peptides": 10,
            "protein_descriptions": ["Description of Protein P12345"],
            "ibaq": 123.45,
            "ribaq": 67.89,
            "intensity": 987.65
        })

    @staticmethod
    def generate_random_protein():
        protein = {
            "protein_accessions": ["P" + str(random.randint(10000, 99999))],
            "sample_accession": "S" + str(random.randint(1, 100)),
            "abundance": round(random.uniform(0.1, 10.0), 2),
            "global_qvalue": round(random.uniform(0.01, 0.1), 2),
            "is_decoy": random.choice([0, 1]),
            "best_id_score": "Score" + str(random.randint(1, 100)) + ": " + str(
                round(random.uniform(0.1, 10.0), 2)),
            "gene_accessions": ["G" + str(random.randint(100, 999)) for _ in range(random.randint(1, 5))],
            "gene_names": ["Gene" + chr(random.randint(65, 90)) for _ in range(random.randint(1, 5))],
            "number_of_peptides": random.randint(1, 50),
            "number_of_psms": random.randint(1, 50),
            "number_of_unique_peptides": random.randint(1, 20),
            "protein_descriptions": ["Description of Protein"],
            "ibaq": round(random.uniform(50.0, 500.0), 2),
            "ribaq": round(random.uniform(10.0, 100.0), 2),
            "intensity": round(random.uniform(100.0, 1000.0), 2)
        }
        return protein

    def test_write_an_example_protein_million(self):

        # Generate a list of random proteins
        num_proteins = 1000000  # You can adjust the number of proteins
        protein_list = [self.generate_random_protein() for _ in range(num_proteins)]

        protein_manager = ProteinHandler()
        protein_manager.parquet_path = 'example_protein.parquet'
        table = protein_manager.create_proteins_table(protein_list)
        protein_manager.write_single_file_parquet(table, write_metadata=True)



    def test_read_parquet_file(self):
        protein_handler = ProteinHandler()
        protein_handler.parquet_path = 'data/example_protein.parquet'
        df_table = protein_handler.read_protein_dataset()
        print(df_table)


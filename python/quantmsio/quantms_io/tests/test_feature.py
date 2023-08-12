from unittest import TestCase

from quantms_io.core.feature import FeatureHandler
import pyarrow as pa
import random


def generate_random_feature():
    """
    Generate a list of features
    """

    fields = FeatureHandler.FEATURE_FIELDS
    num_rows = 1000000
    mock_data = []
    for _ in range(num_rows):
        mock_row = {
            "sequence": "PEPTIDESEQ",
            "protein_accessions": ["P12345", "P67890"],
            "protein_start_positions": [10, 50],
            "protein_end_positions": [30, 70],
            "unique": random.randint(0, 1),
            "modifications": ["Phospho", "Acetyl"],
            "retention_time": random.uniform(40.0, 50.0),
            "charge": random.randint(1, 3),
            "exp_mass_to_charge": random.uniform(2.0, 3.0),
            "calc_mass_to_charge": random.uniform(2.0, 3.0),
            "peptidoform": "PEPTIDE",
            "posterior_error_probability": random.uniform(0.0, 0.1),
            "global_qvalue": random.uniform(0.0, 0.1),
            "is_decoy": random.randint(0, 1),
            "best_id_score": f"Score{random.randint(1, 10)}: {random.uniform(0.0, 1.0)}",
            "intensity": random.uniform(500.0, 1500.0),
            "spectral_count": random.randint(10, 30),
            "sample_accession": "S1",
            "condition": "ConditionA",
            "fraction": "F1",
            "biological_replicate": "BioRep1",
            "fragment_ion": "b2",
            "isotope_label_type": "LabelA",
            "run": "Run1",
            "channel": "ChannelA",
            "id_scores": ["Score1: 0.9", "Score2: 0.6"],
            "consensus_support": random.uniform(0.0, 1.0),
            "reference_file_name": "ref_file.txt",
            "scan_number": "Scan123",
            "mz": [random.uniform(200.0, 800.0) for _ in range(10)],
            "intensity_array": [random.uniform(0.0, 1000.0) for _ in range(10)],
            "num_peaks": random.randint(5, 20),
            "gene_accessions": ["G123", "G456"],
            "gene_names": ["GeneA", "GeneB"]
        }
        mock_data.append(mock_row)
    return mock_data


class TestFeatureHandler(TestCase):
    
    def test_write_an_example_feature_million(self):

        feature_list = generate_random_feature() # Generate a list of features

        feature_manager = FeatureHandler()
        feature_manager.parquet_path = 'example_feature.parquet'
        table = feature_manager.create_feature_table(feature_list)
        feature_manager.write_single_file_parquet(table, write_metadata=True)


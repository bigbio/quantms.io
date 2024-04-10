from unittest import TestCase

from quantms_io.core.tools import generate_features_of_spectrum
from quantms_io.core.tools import map_gene_msgs_to_parquet
from quantms_io.core.tools import map_protein_for_parquet
from quantms_io.core.tools import register_file_to_json


class TestTools(TestCase):

    def test_map_spectrum_message_to_parquet(self):
        mz_foder = __package__ + "/examples/mzml"
        parquet_path = __package__ + "/examples/Tools/test.feature.parquet"
        output_path = __package__ + "/examples/output/Tools/test_fill_mz.feature.parquet"
        label = "feature"
        chunksize = 100000
        partition = "charge"
        generate_features_of_spectrum(parquet_path, mz_foder, output_path, label, chunksize, partition)

    def test_attach_file_to_project_json(self):
        rigister_files = [
            [
                __package__ + path
                for path in [
                    "/examples/output/DDA-plex/MSV000079033-6771dd93-aafd-4530-882c-e98243f53aaf.feature.parquet",
                    "feature_file",
                ]
            ],
            [
                __package__ + path
                for path in [
                    "/examples/output/DDA-plex/MSV000079033-a645fef2-3d54-4166-88a1-519c3be3b817.psm.parquet",
                    "psm_file",
                ]
            ],
        ]
        project_file = __package__ + "/examples/output/DDA-plex/peoject.json"
        for file_item in rigister_files:
            register_file_to_json(
                project_file=project_file, attach_file=file_item[0], category=file_item[1], replace_existing=True
            )

    def test_get_unanimous(self):
        parquet_path = __package__ + "/examples/Tools/test.feature.parquet"
        fasta_path = __package__ + "/examples/fasta/Homo-sapiens.fasta"
        output_path = __package__ + "/examples/output/Tools/test_unanimous.featrue.parquet"
        map_protein_for_parquet(
            parquet_path, fasta_path, output_path, map_parameter="map_protein_name", label="feature"
        )

    def test_generate_gene_msg_to_parquet(self):
        parquet_path = __package__ + "/examples/Tools/test_unanimous.featrue.parquet"
        fasta_path = __package__ + "/examples/fasta/Homo-sapiens.fasta"
        output_path = __package__ + "/examples/output/Tools/test_fill_gene.featrue.parquet"
        label = "feature"
        map_parameter = "map_protein_name"
        species = "human"
        map_gene_msgs_to_parquet(parquet_path, fasta_path, map_parameter, output_path, label, species)

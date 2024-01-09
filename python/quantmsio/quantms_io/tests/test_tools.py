from quantms_io.core.tools import generate_features_of_spectrum
from quantms_io.core.tools import generate_project_report
from quantms_io.core.tools import register_file_to_json

from unittest import TestCase


class TestTools(TestCase):

    def test_map_spectrum_message_to_parquet(self):
        mz_foder = 'data/mzml'
        parquet_path = 'data/PXD014414/PXD014414-943a8f02-0527-4528-b1a3-b96de99ebe75.featrue.parquet'
        output_path = 'data/mzml/PXD014414.parquet'
        label = 'feature'
        chunksize = 100000
        partition = 'charge'
        generate_features_of_spectrum(
            parquet_path,
            mz_foder,
            output_path,
            label,
            chunksize,
            partition
        )
    
    def test_generate_project_report(self):
        project_folder = 'data/PXD014414'
        generate_project_report(project_folder)
        
    def test_attach_file_to_project_json(self):
        rigister_files = [
            ['data/PXD014414/PXD014414-943a8f02-0527-4528-b1a3-b96de99ebe75.featrue.parquet','feature_file'],
            ['data/PXD014414/PXD014414-f4fb88f6-0a45-451d-a8a6-b6d58fb83670.psm.parquet','psm_file'],
            ['data/PXD014414/PXD014414-48fc92e9-ffd2-4335-a6f5-ce7d5278b5ff.differential.tsv','differential_file']
        ]
        project_file = 'data/PXD014414/project.json'
        for file_item in rigister_files:
            register_file_to_json(project_file=project_file,
                                  attach_file=file_item[0],
                                  category=file_item[1],
                                  replace_existing=True
            )
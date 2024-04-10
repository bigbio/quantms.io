import time
from unittest import TestCase

import pandas as pd

from quantms_io.core.diann_convert import DiaNNConvert
from quantms_io.core.diann_convert import find_modification
from quantms_io.core.diann_convert import get_exp_design_dfs
from quantms_io.core.project import create_uuid_filename
from quantms_io.core.tools import generate_start_and_end_from_fasta
from quantms_io.core.tools import load_best_scan_number


class Test(TestCase):
    global test_datas
    test_datas = (
        "/examples/DIA-lfq/diann_report.tsv",
        "/examples/DIA-lfq/mzml",
        "/examples/DIA-lfq/PXD037682.sdrf_openms_design.tsv",
        "/examples/DIA-lfq/PXD037682.sdrf.tsv",
        "/examples/output/DIA-lfq/" + create_uuid_filename("PXD037682", ".feature.parquet"),
        "/examples/output/DIA-lfq/" + create_uuid_filename("PXD037682", ".psm.parquet"),
    )

    def test_all_required_arguments(self):
        st = time.time()

        # Initialize input arguments
        report_path = __package__ + test_datas[0]
        mzml_info_folder = __package__ + test_datas[1]
        design_file = __package__ + test_datas[2]
        sdrf_file = __package__ + test_datas[3]
        feature_output_path = __package__ + test_datas[4]
        psm_output_path = __package__ + test_datas[5]
        qvalue_threshold = 0.01

        DiaNN = DiaNNConvert()

        DiaNN.generate_psm_and_feature_file(
            report_path=report_path,
            design_file=design_file,
            qvalue_threshold=qvalue_threshold,
            feature_output_path=feature_output_path,
            psm_output_path=psm_output_path,
            sdrf_path=sdrf_file,
            mzml_info_folder=mzml_info_folder,
            file_num=10,
            duckdb_max_memory="4GB",
            duckdb_threads=4,
        )
        et = time.time()

        # get the execution time
        elapsed_time = et - st
        print("Execution time:", elapsed_time, "seconds")

    def test_fill_start_and_end(self):
        fasta_path = __package__ + "/examples/fasta/Homo-sapiens.fasta"
        feature_path = __package__ + test_datas[4]
        psm_path = __package__ + test_datas[5]
        feature_output_path = __package__ + "/examples/output/DIA-lfq/PXD010154_fill_start_and_end.feature.parquet"
        psm_output_path = __package__ + "/examples/output/DIA-lfq/PXD010154_fill_start_and_end.psm.parquet"
        generate_start_and_end_from_fasta(
            parquet_path=feature_path, fasta_path=fasta_path, label="feature", output_path=feature_output_path
        )
        generate_start_and_end_from_fasta(
            parquet_path=psm_path, fasta_path=fasta_path, label="psm", output_path=psm_output_path
        )

    def test_fill_best_scan_number(self):
        feature_output_path = __package__ + test_datas[4]
        psm_output_path = __package__ + test_datas[5]
        output_path = __package__ + "/examples/output/DIA-lfq/PXD010154_fill_best_scan_number.feature.parquet"
        load_best_scan_number(
            diann_psm_path=psm_output_path, diann_feature_path=feature_output_path, output_path=output_path
        )

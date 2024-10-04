from .common import datafile
from ddt import data
from ddt import ddt
from unittest import TestCase
from quantmsio.operate.query import Query

@ddt
class TestHandler(TestCase):
    global test_datas
    psm_path = datafile("parquet/psm.parquet")
    feature_path = datafile("parquet/feature.parquet")
    fasta = datafile("fasta/Homo-sapiens.fasta")
    test_datas = [
        (
            "/examples/mzml/test.parquet",
            "/examples/mzml/mzml",
        ),
    ]

    def test_iter_samples(self):
        q = Query(TestHandler.feature_path)
        for _ in q.iter_samples():
            print("ok")
    
    def test_iter_chunk(self):
        q = Query(TestHandler.feature_path)
        for _ in q.iter_chunk():
            print("ok")
    
    def test_iter_file(self):
        q = Query(TestHandler.feature_path)
        for _ in q.iter_file():
            print("ok")
    
    @data(*test_datas)
    def test_inject_spectrum_msg(self,test_data):
        parquet_file = datafile(test_data[0])
        mz_folder = datafile(test_data[1])
        q = Query(parquet_file)
        df = q.get_report_from_database(["test"])
        df["scan_number"] = "1"
        q.inject_spectrum_msg(df, mz_folder)

    def test_inject_position_msg(self):
        q = Query(TestHandler.feature_path)
        df = q.get_report_from_database(["20180914_QE8_nLC0_BDA_SA_DIA_Keratinocytes_NN002"])
        protein_dict = q.get_protein_dict(TestHandler.fasta)
        q.inject_position_msg(df,protein_dict)
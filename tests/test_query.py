from .common import datafile
from unittest import TestCase
from quantmsio.operate.query import Query


class TestHandler(TestCase):
    feature_path = datafile("parquet/feature.parquet")
    fasta = datafile("fasta/Homo-sapiens.fasta")

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

    def test_inject_position_msg(self):
        q = Query(TestHandler.feature_path)
        df = q.get_report_from_database(["20180914_QE8_nLC0_BDA_SA_DIA_Keratinocytes_NN002"])
        protein_dict = q.get_protein_dict(TestHandler.fasta)
        q.inject_position_msg(df, protein_dict)

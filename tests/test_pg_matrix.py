from .common import datafile
from ddt import ddt
from unittest import TestCase
from quantmsio.core.pg_matrix import PgMatrix


@ddt
class TestHandler(TestCase):
    feature_path = datafile("parquet/feature.parquet")

    def test_iter_samples(self):
        q = PgMatrix(TestHandler.feature_path)
        for _ in q.generate_pg_matrix():
            print("ok")

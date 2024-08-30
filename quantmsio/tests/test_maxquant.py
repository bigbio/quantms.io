from unittest import TestCase

from ddt import data
from ddt import ddt
from quantmsio.core.maxquant_convert import MaxquantConvert
from quantmsio.core.maxquant_convert import get_mods, get_mod_map
from quantmsio.tests.common import datafile

@ddt
class TestMaxquantHandler(TestCase):
    global test_datas
    test_datas = [
        (
            "MaxQuant-lfq/example_sdrf.tsv",
            "MaxQuant-lfq/evidence.txt",
        ),
    ]

    @data(*test_datas)
    def test_maxquant_convert(self, test_data):
        C  = MaxquantConvert()
        C._modifications = get_mods(datafile(test_data[0]))
        mods_map = get_mod_map(datafile(test_data[0]))
        for df in C.iter_batch(datafile(test_data[1]),mods_map,50000):
             C.merge_sdrf(df,datafile(test_data[0]))


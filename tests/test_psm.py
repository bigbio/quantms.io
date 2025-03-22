from .common import datafile
from quantmsio.core.psm import Psm


def test_convert_mztab_to_feature():
    mztab_path = datafile("DDA-lfq/PXD040438.mzTab")
    psm = Psm(mztab_path)
    for _ in psm.generate_report():
        print("ok")

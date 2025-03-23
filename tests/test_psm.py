from pathlib import Path

from quantmsio.core.psm import Psm

TEST_DATA_ROOT = Path(__file__).parent / "examples"


def test_convert_mztab_to_feature():
    mztab_path = TEST_DATA_ROOT / "DDA-lfq/PXD040438.mzTab"
    psm = Psm(mztab_path)
    for _ in psm.generate_report():
        print("ok")

from pathlib import Path
from quantmsio.operate.query import Query

TEST_DATA_ROOT = Path(__file__).parent / "examples"
feature_path = TEST_DATA_ROOT / "parquet/feature.parquet"
fasta = TEST_DATA_ROOT / "fasta/Homo-sapiens.fasta"


def test_iter_samples():
    q = Query(feature_path)
    for _ in q.iter_samples():
        print("ok")


def test_iter_chunk():
    q = Query(feature_path)
    for _ in q.iter_chunk():
        print("ok")


def test_iter_file():
    q = Query(feature_path)
    for _ in q.iter_file():
        print("ok")


def test_inject_position_msg():
    q = Query(feature_path)
    df = q.get_report_from_database(
        ["20180914_QE8_nLC0_BDA_SA_DIA_Keratinocytes_NN002"]
    )
    protein_dict = q.get_protein_dict(fasta)
    q.inject_position_msg(df, protein_dict)

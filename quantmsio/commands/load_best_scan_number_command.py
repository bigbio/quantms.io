import click

from quantmsio.core.tools import load_best_scan_number


@click.command(
    "inject-bset-psm-scan-number",
    short_help="inject bset_psm_scan_number to feature",
)
@click.option("--diann_psm_path", help="diann psm parquet file path", required=True)
@click.option("--diann_feature_path", help="diann feature parquet file path", required=True)
@click.option("--output_path", help="save path", required=True)
def inject_bset_psm_scan_number(diann_psm_path: str, diann_feature_path: str, output_path: str):
    """
    Register the file with project.json
    :param diann_psm_path: diann psm parquet file path
    :param diann_feature_path: diann feature parquet file path
    :param output_path: save path
    :return: none
    """
    load_best_scan_number(
        diann_psm_path=diann_psm_path,
        diann_feature_path=diann_feature_path,
        output_path=output_path,
    )

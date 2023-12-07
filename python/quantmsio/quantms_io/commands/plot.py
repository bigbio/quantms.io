import click
from quantms_io.core.tools import plot_peptides_of_lfq_condition,plot_distribution_of_ibaq
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


# Command Group
@click.group(name='plot', context_settings=CONTEXT_SETTINGS)
def plot():
    """Tool related commands"""
    pass

@plot.command("plot-peptides", short_help="plot peptides of condition in lfq",)
@click.option("--psm_parquet_path", help="psm parquet path in lfq", required=True)
@click.option("--sdrf_path", help="sdrf path", required=True)
@click.option("--save_path", help="img save path [xxx.svg]", required=True)
@click.pass_context
def plot_peptides(ctx, psm_parquet_path: str, sdrf_path: str, save_path: str):
    """
    convert mztab psm section to a parquet file. The parquet file will contain the features and the metadata.
    :param psm_parquet_path: psm parquet path in lfq
    :param sdrf_path: sdrf path
    :param save_path: img save path [xxx.png]
    :return: none
    """
    plot_peptides_of_lfq_condition(psm_parquet_path=psm_parquet_path,sdrf_path=sdrf_path,save_path=save_path)

@plot.command("plot-ibaq-distribution", short_help="plot ibaq distribution of expression")
@click.option('--ibaq_path',  help='ibaq file path',required=True)
@click.option('--save_path',  help='img save path [xxx.svg]',required=True)
@click.pass_context
def plot_ibaq_distribution(ctx,ibaq_path,save_path):
    """
    plot ibaq distribution of expression
    :param ibaq_path: ibaq file path
    :param save_path: img save path [xxx.png]
    :return: none
    """

    plot_distribution_of_ibaq(ibaq_path,save_path)
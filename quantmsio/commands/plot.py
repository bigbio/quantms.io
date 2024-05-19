import click

from quantmsio.core.plots import plot_distribution_of_ibaq
from quantmsio.core.plots import plot_intensity_box_of_samples
from quantmsio.core.plots import plot_intensity_distribution_of_samples
from quantmsio.core.plots import plot_peptide_distribution_of_protein
from quantmsio.core.plots import plot_peptides_of_lfq_condition

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


# Command Group
@click.group(name="plot", context_settings=CONTEXT_SETTINGS)
def plot():
    """Tool related commands"""
    pass


@plot.command(
    "plot-psm-peptides",
    short_help="plot peptides of condition in lfq",
)
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
    plot_peptides_of_lfq_condition(psm_parquet_path=psm_parquet_path, sdrf_path=sdrf_path, save_path=save_path)


@plot.command("plot-ibaq-distribution", short_help="plot ibaq distribution of expression")
@click.option("--ibaq_path", help="ibaq file path", required=True)
@click.option("--save_path", help="img save path [xxx.svg]", required=True)
@click.option("--select_column", help="Selected column in Ibaq File", required=False)
@click.pass_context
def plot_ibaq_distribution(ctx, ibaq_path: str, save_path: str, select_column: str) -> None:
    """
    plot ibaq distribution of expression
    :param ibaq_path: ibaq file path
    :param save_path: img save path [xxx.png]
    :return: none
    """

    plot_distribution_of_ibaq(ibaq_path, save_path, select_column)


@plot.command(
    "plot-kde-intensity-distribution",
    short_help="plot kde intensity distribution of feature",
)
@click.option("--feature_path", help="feature file path", required=True)
@click.option("--save_path", help="img save path [xxx.svg]", required=True)
@click.option("--num_samples", help="The number of samples plotted", default=10)
@click.pass_context
def plot_kde_intensity_distribution(feature_path: str, save_path: str, num_samples: int):
    """
    plot ibaq distribution of expression
    :param feature_path: feature file path
    :param save_path: img save path [xxx.png]
    :param num_samples: The number of samples plotted
    :return: none
    """

    plot_intensity_distribution_of_samples(feature_path, save_path, num_samples)


@plot.command(
    "plot-bar-peptide-distribution",
    short_help="plot bar peptide distribution of feature",
)
@click.option("--feature_path", help="feature file path", required=True)
@click.option("--save_path", help="img save path [xxx.svg]", required=True)
@click.option("--num_samples", help="The number of samples plotted", default=20)
@click.pass_context
def plot_bar_peptide_distribution(feature_path: str, save_path: str, num_samples: int):
    """
    plot ibaq distribution of expression
    :param feature_path: feature file path
    :param save_path: img save path [xxx.png]
    :param num_samples: The number of samples plotted
    :return: none
    """

    plot_peptide_distribution_of_protein(feature_path, save_path, num_samples)


@plot.command(
    "plot-box-intensity-distribution",
    short_help="plot box intensity distribution of feature",
)
@click.option("--feature_path", help="feature file path", required=True)
@click.option("--save_path", help="img save path [xxx.svg]", required=True)
@click.option("--num_samples", help="The number of samples plotted", default=10)
@click.pass_context
def plot_box_intensity_distribution(feature_path: str, save_path: str, num_samples: int):
    """
    plot ibaq distribution of expression
    :param feature_path: feature file path
    :param save_path: img save path [xxx.png]
    :param num_samples: The number of samples plotted
    :return: none
    """
    plot_intensity_box_of_samples(feature_path, save_path, num_samples)

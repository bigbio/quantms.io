import sys

import click

from quantmsio.core.statistics import IbaqStatistics
from quantmsio.core.statistics import ParquetStatistics

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


# Command Group
@click.group(name="statistics", context_settings=CONTEXT_SETTINGS)
def statistics():
    """Tool related commands"""
    pass


@statistics.command(
    "project-ae-statistics",
    short_help="Statistics about a particular project",
)
@click.option("--absolute_path", help="absolute path", required=True)
@click.option("--parquet_path", help="psm parquet path in lfq", required=True)
@click.option(
    "--save_path",
    help="file with the statistics (e.g. statistics.csv), if not provided," " will print to stdout",
)
@click.pass_context
def feature_file_statistics(ctx, absolute_path: str, parquet_path: str, save_path: str):
    """
    Statistics of a feature file
    :param parquet_path: feature parquet path
    :param save_path: file with the statistics (e.g. statistics.csv), if not provided, will print to stdout
    :return: none
    """
    feature_statistics = ParquetStatistics(parquet_path)
    absolute_stats = IbaqStatistics(ibaq_path=absolute_path)

    def write_stats(file, stats: ParquetStatistics):
        file.write("Number of proteins: {}\n".format(stats.get_number_of_proteins()))
        file.write("Number of peptides: {}\n".format(stats.get_number_of_peptides()))
        file.write("Number of samples: {}\n".format(stats.get_number_of_samples()))
        file.write("Number of peptidoforms: {}\n".format(stats.get_number_of_peptidoforms()))
        file.write("Number of msruns: {}\n".format(stats.get_number_msruns()))

    def write_absolute_stats(file, stats: IbaqStatistics):
        file.write("Ibaq Number of proteins: {}\n".format(stats.get_number_of_proteins()))
        file.write("Ibaq Number of samples: {}\n".format(stats.get_number_of_samples()))

    if save_path:
        # Open save file and write stats
        with open(save_path, "w") as f:
            write_stats(f, feature_statistics)
            write_absolute_stats(f, absolute_stats)
    else:
        # Print stats to stdout
        write_stats(sys.stdout, feature_statistics)
        write_absolute_stats(sys.stdout, absolute_stats)


@statistics.command(
    "parquet-psm-statistics",
    short_help="Statistics about a particular psm parquet file",
)
@click.option("--parquet_path", help="psm parquet path in lfq", required=True)
@click.option(
    "--save_path",
    help="file with the statistics (e.g. statistics.csv), if not provided," " will print to stdout",
)
@click.pass_context
def parquet_psm_statistics(ctx, parquet_path: str, save_path: str):
    """
    Statistics of a psm parquet file
    :param parquet_path: psm parquet path
    :param save_path: file with the statistics (e.g. statistics.csv), if not provided, will print to stdout
    :return: none
    """

    def write_stats(file, stats: ParquetStatistics):
        file.write("Number of proteins: {}\n".format(stats.get_number_of_proteins()))
        file.write("Number of peptides: {}\n".format(stats.get_number_of_peptides()))
        file.write("Number of peptidoforms: {}\n".format(stats.get_number_of_peptidoforms()))
        file.write("Number of psms: {}\n".format(stats.get_number_of_psms()))
        file.write("Number of msruns: {}\n".format(stats.get_number_msruns()))

    feature_statistics = ParquetStatistics(parquet_path)
    if save_path:
        # Open save file and write stats
        with open(save_path, "w") as f:
            write_stats(f, feature_statistics)
    else:
        # Print stats to stdout
        write_stats(sys.stdout, feature_statistics)

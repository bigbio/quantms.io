import sys

import click

from quantms_io.core.statistics import ParquetStatistics

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


# Command Group
@click.group(name='statistics', context_settings=CONTEXT_SETTINGS)
def plot():
    """Tool related commands"""
    pass

@plot.command("feature-file-statistics", short_help="Statistics of a feature file",)
@click.option("--parquet_path", help="psm parquet path in lfq", required=True)
@click.option("--save_path", help="file with the statistics (e.g. statistics.csv), if not provided,"
                                  " will print to stdout")
@click.pass_context
def feature_file_statistics(ctx, parquet_path: str, save_path: str):
    """
    Statistics of a feature file
    :param parquet_path: feature parquet path
    :param save_path: file with the statistics (e.g. statistics.csv), if not provided, will print to stdout
    :return: none
    """
    feature_statistics = ParquetStatistics(parquet_path)

    def write_stats(file, stats):
        file.write("Number of proteins: {}\n".format(stats.get_number_of_proteins()))
        file.write("Number of peptides: {}\n".format(stats.get_number_of_peptides()))
        file.write("Number of samples: {}\n".format(stats.get_number_of_samples()))
        file.write("Number of peptidoforms: {}\n".format(stats.get_number_of_peptidoforms()))
        file.write("Number of msruns: {}\n".format(stats.get_number_msruns()))

    if save_path:
        # Open save file and write stats
        with open(save_path, 'w') as f:
            write_stats(f, feature_statistics)
    else:
        # Print stats to stdout
        write_stats(sys.stdout, feature_statistics)

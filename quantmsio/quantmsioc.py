"""
Commandline interface for quantmsio package allows generating the quantms.io file format from different sources and
 steps. The quantms.io specification is available in the docs folder of this repository.
"""

import logging
import click

from quantmsio import __version__ as __version__
from quantmsio.commands.project_command import generate_pride_project_json
from quantmsio.commands.feature_command import convert_feature_file
from quantmsio.commands.psm_command import convert_psm_file, compare_set_of_psms
from quantmsio.commands.diann_command import (
    diann_convert_to_parquet,
    diann_pg_convert_to_parquet,
)
from quantmsio.commands.ae_command import convert_ibaq_absolute
from quantmsio.commands.de_command import convert_msstats_differential
from quantmsio.commands.attach_file_command import attach_file_to_json
from quantmsio.commands.generate_spectra_message_command import (
    map_spectrum_message_to_parquet,
)
from quantmsio.commands.generate_gene_message_command import map_gene_message_to_parquet
from quantmsio.commands.plot_command import plot
from quantmsio.commands.statistic_command import statistics
from quantmsio.commands.maxquant_command import (
    convert_maxquant_psm,
    convert_maxquant_feature,
)
from quantmsio.commands.ibaq_command import convert_ibaq_file
from quantmsio.commands.fragpipe_command import convert_fragpipe_psm
from quantmsio.commands.map_latest_uniport_command import map_latest_uniport
from quantmsio.commands.anndata_command import merge_ae_files

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.version_option(
    version=__version__, package_name="quantmsio", message="%(package)s %(version)s"
)
@click.group(context_settings=CONTEXT_SETTINGS)
def cli() -> None:
    """
    This is the main tool that gives access to all commands to convert SDRF files into pipeline-specific configuration files
    """
    logging.basicConfig(
        level=logging.INFO,
        datefmt="%H:%M:%S",
        format="[%(asctime)s] %(levelname).1s | %(name)s | %(message)s",
    )


cli.add_command(generate_pride_project_json)
cli.add_command(convert_feature_file)
cli.add_command(convert_psm_file)
cli.add_command(compare_set_of_psms)
cli.add_command(diann_convert_to_parquet)
cli.add_command(diann_pg_convert_to_parquet)
cli.add_command(convert_ibaq_absolute)
cli.add_command(convert_msstats_differential)
cli.add_command(attach_file_to_json)
cli.add_command(map_spectrum_message_to_parquet)
cli.add_command(map_gene_message_to_parquet)
cli.add_command(plot)
cli.add_command(statistics)
cli.add_command(convert_maxquant_psm)
cli.add_command(convert_maxquant_feature)
cli.add_command(convert_ibaq_file)
cli.add_command(map_latest_uniport)
cli.add_command(convert_fragpipe_psm)
cli.add_command(merge_ae_files)


def quantms_io_main() -> None:
    """
    Main function to run the quantmsio command line interface
    :return: none
    """
    cli()


if __name__ == "__main__":
    quantms_io_main()

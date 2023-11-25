import click

from quantms_io.core.de import DifferentialExpressionHandler
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    This is the main tool that gives access to all commands.
    """

@click.command("convert-de", short_help="Convert a MSstats differential file into a quantms.io file format")
@click.option("--msstats_file", help="MSstats differential file", required=True)
@click.option(
    "--sdrf_file",
    help="the SDRF file needed to extract some of the metadata",
    required=True,
)
@click.option(
    "--project_file",
    help="quantms.io project file",
    required=False,
)
@click.option(
    "--fdr_threshold",
    help="FDR threshold to use to filter the results",
    required=False,
    default="0.05",
)
@click.option(
    "--output_folder", help="Folder to generate the df expression file.", required=True
)
@click.option("--output_prefix_file", help="Prefix of the df expression file", required=False)
@click.option(
    "--delete_existing", help="Delete existing files in the output folder", is_flag=True
)
@click.pass_context
def convert_msstats_differential(
    ctx,
    msstats_file: str,
    sdrf_file: str,
    project_file:str,
    fdr_threshold: float,
    output_folder: str,
    output_prefix_file: str,
    delete_existing: bool=True,
):
    """
    Convert a MSstats differential file into a quantms.io file format. The file definition is available in the docs
    https://github.com/bigbio/quantms.io/blob/main/docs/DE.md.
    :param msstats_file: MSstats differential file
    :param sdrf_file: the SDRF file needed to extract some of the metadata
    :param project_file: quantms.io project file
    :param output_folder: Folder to generate the df expression file.
    :param output_prefix_file: Prefix of the df expression file
    :param delete_existing: Delete existing files in the output folder
    :param fdr_threshold: FDR threshold to use to filter the results
    :return: none
    """

    if msstats_file is None or sdrf_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    de_handler = DifferentialExpressionHandler()
    if project_file:
        de_handler.load_project_file(project_file)
    de_handler.load_msstats_file(msstats_file)
    de_handler.load_sdrf_file(sdrf_file)
    de_handler.set_fdr_threshold(fdr_threshold=fdr_threshold)
    de_handler.convert_msstats_to_quantms(
        output_folder=output_folder,
        output_file_prefix=output_prefix_file,
        delete_existing=delete_existing,
    )

cli.add_command(convert_msstats_differential)
if __name__ == '__main__':
    cli()
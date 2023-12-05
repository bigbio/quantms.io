from quantms_io.core.ae import AbsoluteExpressionHander
import click

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    This is the main tool that gives access to all commands.
    """


@click.command(
    "convert-ae",
    short_help="Convert a ibaq_absolute file into a quantms.io file " "format",
)
@click.option(
    "--ibaq_file",
    help="the ibaq file path",
    required=True,
)
@click.option(
    "--sdrf_file",
    help="the sdrf_file path",
    required=True,
)
@click.option(
    "--project_file",
    help="quantms.io project file",
    required=False,
)
@click.option(
    "--output_folder", help="Folder to generate the df expression file.", required=True
)
@click.option("--output_prefix_file", help="Prefix of the df expression file", required=False)
@click.option(
    "--delete_existing", help="Delete existing files in the output folder", is_flag=True
)
@click.pass_context
def convert_ibaq_absolute(ctx,
                          ibaq_file: str,
                          sdrf_file: str,
                          project_file: str,
                          output_folder: str,
                          output_prefix_file: str,
                          delete_existing: bool = True,
                          ):
    """
    Convert a IBAQ absolute file into a quantms.io file format. The file definition is available in the docs
    https://github.com/bigbio/quantms.io/blob/main/docs/AE.md.
    :param ibaq_file: IBAQ file
    :param sdrf_file: sdrf file
    :param project_file: quantms.io project file
    :param output_folder: Folder to generate the df expression file.
    :param output_prefix_file: Prefix of the df expression file
    :param delete_existing: Delete existing files in the output folder
    :return: none
    """
    ae_handler = AbsoluteExpressionHander()
    if project_file:
        ae_handler.load_project_file(project_file)
    ae_handler.load_ibaq_file(ibaq_file)
    ae_handler.load_sdrf_file(sdrf_file)
    ae_handler.convert_ibaq_to_quantms(
        output_folder=output_folder,
        output_file_prefix=output_prefix_file,
        delete_existing=delete_existing,
    )


cli.add_command(convert_ibaq_absolute)

if __name__ == '__main__':
    cli()

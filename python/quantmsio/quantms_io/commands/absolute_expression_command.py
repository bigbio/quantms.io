from quantms_io.core.ae import AbsoluteExpressionHander
import click
@click.command(
    "convert-ibaq-absolute",
    short_help="Convert a ibaq_absolute file into a quantms.io file " "format",
)
@click.option(
    "--ibaq_file",
    help="the ibaq file path",
    required=True,
)
@click.option(
    "--output_folder", help="Folder to generate the df expression file.", required=True
)
@click.option("--output_file", help="Prefix of the df expression file", required=False)
@click.option(
    "--generate_project",
    help="Generate project.json for pride project, Otherwise, False",
    required=False,
    is_flag=True,
)
@click.option(
    "--delete_existing", help="Delete existing files in the output folder", is_flag=True
)
@click.option("--project_file", help="quantms.io project file", required=False)
def convert_ibaq_absolute(
    ibaq_file:str,
    output_folder: str,
    output_file: str,
    generate_project: bool = True,
    delete_existing: bool = True,
    project_file:str = None,
    ):
    """
    Convert a IBAQ absolute file into a quantms.io file format. The file definition is available in the docs
    https://github.com/bigbio/quantms.io/blob/main/docs/AE.md.
    :param ibaq_file: IBAQ file
    :param project_file: quantms.io project file
    :param output_folder: Folder to generate the df expression file.
    :param output_file: Prefix of the df expression file
    :param generate_project: "Generate project.json for pride project, Otherwise, False"
    :param delete_existing: Delete existing files in the output folder
    :return: none
    """
    de_handler = AbsoluteExpressionHander()
    if generate_project:
        de_handler.load_project_file(project_file)
    de_handler.load_ibaq_file(ibaq_file)
    de_handler.convert_ibaq_to_quantms(
        output_folder=output_folder,
        output_file_prefix=output_file,
        delete_existing=delete_existing,
    )
    if generate_project:
        de_handler.update_project_file(project_file)
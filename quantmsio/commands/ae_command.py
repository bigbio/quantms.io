import click

from quantmsio.core.ae import AbsoluteExpressionHander
from quantmsio.utils.file_utils import extract_protein_list


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
    "--protein_file",
    help="Protein file that meets specific requirements",
    required=False,
)
@click.option(
    "--project_file",
    help="quantms.io project file",
    required=False,
)
@click.option(
    "--output_folder", help="Folder to generate the df expression file.", required=True
)
@click.option(
    "--output_prefix_file", help="Prefix of the df expression file", required=False
)
@click.option(
    "--delete_existing", help="Delete existing files in the output folder", is_flag=True
)
def convert_ibaq_absolute(
    ibaq_file: str,
    sdrf_file: str,
    project_file: str,
    protein_file: str,
    output_folder: str,
    output_prefix_file: str,
    delete_existing: bool = True,
):
    protein_list = extract_protein_list(protein_file) if protein_file else None
    protein_str = "|".join(protein_list) if protein_list else None
    ae_handler = AbsoluteExpressionHander()
    if project_file:
        ae_handler.load_project_file(project_file)
    ae_handler.load_ibaq_file(ibaq_file, protein_str)
    ae_handler.load_sdrf_file(sdrf_file)
    ae_handler.convert_ibaq_to_quantms(
        output_folder=output_folder,
        output_file_prefix=output_prefix_file,
        delete_existing=delete_existing,
    )

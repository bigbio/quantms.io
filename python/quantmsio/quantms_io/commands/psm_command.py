import click

from quantms_io.core.psm import PSMHandler
from quantms_io.core.project import create_uuid_filename,check_directory
from quantms_io.core.tools import plot_peptidoform_charge_venn, plot_sequence_venn


@click.command(
    "convert-psm-file",
    short_help="Convert psm from mzTab to parquet file in quantms io",
)
@click.option(
    "--mztab_file",
    help="the mzTab file, this will be used to extract the protein information",
    required=True,
)
@click.option(
    "--output_folder",
    help="Folder where the parquet file will be generated",
    required=True,
)
@click.option(
    "--generate_project",
    help="Generate project.json for pride project, Otherwise, False",
    required=False,
    is_flag=True,
)
@click.option(
    "--output_prefix_file",
              help="Prefix of the parquet file needed to generate the file name",
              required=False)
@click.option(
    "--verbose",
              help="Output debug information.",
    default=False, is_flag=True)
def convert_psm_file(mztab_file: str, output_folder: str,generate_project:bool = True, output_prefix_file: str=None, verbose: bool = False):
    """
    Convert mztab psm section to a parquet file. The parquet file will contain the features and the metadata.
    :param mztab_file: the mzTab file, this will be used to extract the protein information
    :param output_folder: Folder where the Json file will be generated
    :param generate_project: "Generate project.json for pride project, Otherwise, False"
    :param output_prefix_file: Prefix of the Json file needed to generate the file name
    :param verbose: Output debug information.
    :return: none
    """

    if mztab_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    if generate_project:
        Project = check_directory(output_folder)
        project_accession = Project.project.project_info["project_accession"]

    if not output_prefix_file:
        if generate_project:
            output_prefix_file = project_accession
        else:
            output_prefix_file = ''

    psm_manager = PSMHandler()
    psm_manager.parquet_path = output_folder + "/" + create_uuid_filename(output_prefix_file,'.psm.parquet')
    psm_manager.convert_mztab_to_psm(
        mztab_path=mztab_file, output_folder=output_folder, parquet_path=psm_manager.parquet_path, verbose=verbose, generate_project = generate_project
    )


@click.command(
    "compare-set-of-psms", short_help="plot venn for a set of Psms parquet"
)
@click.option('--parquets', type=str, help='List of psm parquet path', multiple=True)
@click.option('--tags', type=str, help='List of parquet label', multiple=True)
def compare_set_of_psms(parquets, tags):
    """
    Compare a set of psm parquet files
    :param parquets: a set of psm parquet path
    :param tags: a set of psm label
    """
    if len(parquets) != len(tags):
        raise click.UsageError("Please provide same length of parquet_list and label_list")

    plot_peptidoform_charge_venn(parquets, tags)
    plot_sequence_venn(parquets, tags)

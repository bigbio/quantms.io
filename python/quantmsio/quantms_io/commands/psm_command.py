import click

from quantms_io.core.psm import PSMHandler
from quantms_io.core.project import create_uuid_filename,check_directory

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
@click.option("--output_prefix_file", help="Prefix of the parquet file needed to generate the file name", required=False)
def convert_psm_file(mztab_file: str, output_folder: str,output_prefix_file: str=None):
    """
    Convert mztab psm section to a parquet file. The parquet file will contain the features and the metadata.
    :param mztab_file: the mzTab file, this will be used to extract the protein information
    :param output_folder: Folder where the Json file will be generated
    :param output_prefix_file: Prefix of the Json file needed to generate the file name
    :return: none
    """

    if mztab_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    Project = check_directory(output_folder)
    project_accession = Project.project.project_info["project_accession"]
    if not output_prefix_file:
        output_prefix_file = project_accession

    psm_manager = PSMHandler()
    psm_manager.parquet_path = output_folder + "/" + create_uuid_filename(output_prefix_file,'.psm.parquet')
    psm_manager.convert_mztab_to_psm(
        mztab_path=mztab_file, parquet_path=psm_manager.parquet_path
    )

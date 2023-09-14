import click

from quantms_io.core.project import ProjectHandler
from quantms_io.core.project import check_directory

@click.command(
    "generate-pride-project-json",
    short_help="Generate a json project file from original pride accession",
)
@click.option(
    "--project_accession",
    help="Project accession to generate the json file",
    required=True,
)
@click.option(
    "--sdrf",
    help="the SDRF file needed to extract some of the metadata",
    required=True,
)
@click.option(
    "--quantms_version",
    help="Quantms.io version to use as output in the project file",
    required=True,
)
@click.option(
    "--output_folder",
    help="Folder where the Json file will be generated",
    required=True,
)

#@click.option(
#    "--output_file",
#    help="Prefix of the Json file needed to generate the file name",
#    required=False,
#)
#@click.option(
#    "--delete_existing", help="Delete existing files in the output folder", is_flag=True
#)
def generate_pride_project_json(
    project_accession: str,
    sdrf: str,
    quantms_version: str,
    output_folder: str,
    #output_file: str,
    #delete_existing: bool,
):
    """
    Generate a json project file from the original PX accession and SDRF file. The Json file definition is available in the docs
    folder of this repository https://github.com/bigbio/quantms.io/blob/main/docs/PROJECT.md. This command will generate
    the Json and attach to the project files the SDRF file provided as input.

    :param project_accession: Project accession to generate the json file
    :param sdrf_file: the SDRF file needed to extract some of the metadata
    :param quantms_version: Quantms.io version to use as output in the project file
    :param output_folder: Folder where the Json file will be generated
    :param output_file: Prefix of the Json file needed to generate the file name
    :param delete_existing: Delete existing files in the output folder
    :return: none
    """

    if (
        project_accession is None
        or sdrf is None
        or quantms_version is None
        or output_folder is None
    ):
        raise click.UsageError("Please provide all the required parameters")
    
    

    #if output_file is None:
        #output_file = project_accession

    project_handler = check_directory(output_folder,project_accession)

    # Populate the project handler with the metadata from Pride Archive and the SDRF file
    project_handler.populate_from_pride_archive()
    project_handler.populate_from_sdrf(sdrf)
    project_handler.add_quantms_version(quantms_version=quantms_version)
    project_path = output_folder + '/' + 'project.json'
    project_handler.save_updated_project_info(output_file_name=project_path)
    '''
    project_handler.add_sdrf_file(
        sdrf_file_path=sdrf_file,
        output_folder=output_folder,
        delete_existing=delete_existing,
    )
    project_handler.save_project_info(
        output_prefix_file=output_file,
        output_folder=output_folder,
        delete_existing=delete_existing,
    )
    '''

import click

from quantmsio.core.project import ProjectHandler


@click.command(
    "attach-file",
    short_help="Register the file to project.json.",
)
@click.option("--project_file", help="the project.json file", required=True)
@click.option("--attach_file", help="The path of the file that will be registered", required=True)
@click.option(
    "--category",
    type=click.Choice(
        ["sdrf_file", "feature_file", "psm_file", "differential_file", "absolute_file"],
        case_sensitive=False,
    ),
    help="The type of file that will be registered",
    required=True,
)
@click.option("--is_folder", help="A boolean value that indicates if the file is a folder or not", is_flag=True)
@click.option(
    "--partitions",
    help="The field used for splitting files, multiple fields are separated by ,",
    required=False,
)
@click.option("--replace_existing", help="Whether to delete old files", is_flag=True)
def attach_file_to_json(project_file, attach_file, category, is_folder, partitions, replace_existing):
    """
    Register the file with project.json
    :param path_name: The name of the file or folder
    :param file_category: quantms file category(e.g."protein_file","peptide_file","psm_file","differential_file",etc.)
    :param is_folder: A boolean value that indicates if the file is a folder or not.
    :partitions: The fields that are used to partition the data in the file. This is used to optimize the data retrieval and filtering of the data. This field is optional.
    :param replace_existing: Whether to delete old files
    """
    if partitions:
        partitions = partitions.split(",")
    register = ProjectHandler(project_json_file=project_file)
    register.register_file(
        attach_file,
        category,
        is_folder=is_folder,
        partition_fields=partitions,
        replace_existing=replace_existing,
    )
    register.save_updated_project_info(output_file_name=project_file)

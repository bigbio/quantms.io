import click
from quantmsio.core.project import create_uuid_filename
from quantmsio.core.pg_matrix import PgMatrix


@click.command(
    "convert-pg-matrix",
    short_help="Convert psm from mzTab to parquet file in quantms io",
)
@click.option(
    "--feature_file",
    help="feature file",
    required=True,
)
@click.option(
    "--output_folder",
    help="Folder where the parquet file will be generated",
    required=True,
)
@click.option(
    "--output_prefix_file",
    help="Prefix of the parquet file needed to generate the file name",
    required=False,
)
def convert_pg_matrix(
    feature_file: str,
    output_folder: str,
    output_prefix_file: str,
):
    """
    :param feature_file: feature file
    :param output_folder: Folder where the Json file will be generated
    :param output_prefix_file: Prefix of the Json file needed to generate the file name
    """

    if feature_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    if not output_prefix_file:
        output_prefix_file = ""

    pg_manager = PgMatrix(feature_path=feature_file)
    output_path = output_folder + "/" + create_uuid_filename(output_prefix_file, ".pg.parquet")
    pg_manager.write_to_file(output_path=output_path)

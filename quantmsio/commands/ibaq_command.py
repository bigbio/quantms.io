import click
from quantmsio.core.project import create_uuid_filename
from quantmsio.operate.tools import write_ibaq_feature


@click.command(
    "convert-ibaq",
    short_help="Convert psm from mzTab to parquet file in quantms io",
)
@click.option(
    "--feature_file",
    help="feature file",
    required=True,
)
@click.option(
    "--sdrf_file",
    help="the SDRF file needed to extract some of the metadata",
    required=True,
)
@click.option(
    "--output_folder",
    help="Folder where the Json file will be generated",
    required=True,
)
@click.option(
    "--output_prefix_file",
    help="Prefix of the parquet file needed to generate the file name",
    required=False,
)
def convert_ibaq_file(
    feature_file: str,
    sdrf_file: str,
    output_folder: str,
    output_prefix_file: str,
):
    """
    :param feature_file: feature file
    :param sdrf_file: the SDRF file needed to extract some of the metadata
    :param output_folder: Folder where the Json file will be generated
    :param output_prefix_file: Prefix of the Json file needed to generate the file name
    """

    if feature_file is None or sdrf_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    if not output_prefix_file:
        output_prefix_file = ""

    output_path = (
        output_folder + "/" + create_uuid_filename(output_prefix_file, ".ibaq.parquet")
    )
    write_ibaq_feature(sdrf_file, feature_file, output_path)

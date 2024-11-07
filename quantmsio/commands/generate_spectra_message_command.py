import click
from quantmsio.operate.tools import generate_features_of_spectrum


@click.command(
    "map-spectrum-message-to-parquet",
    short_help="According mzMl to map the spectrum message to parquet",
)
@click.option("--parquet_path", help="Psm or feature parquet path")
@click.option("--mzml_directory", help="mzml file folder")
@click.option(
    "--output_folder",
    help="Folder where the Json file will be generated",
    required=True,
)
@click.option(
    "--file_num",
    help="The number of rows of parquet read using pandas streaming",
    default=10,
)
@click.option(
    "--partitions",
    help="The field used for splitting files, multiple fields are separated by ,",
    required=False,
)
def map_spectrum_message_to_parquet(
    parquet_path: str,
    mzml_directory: str,
    output_folder: str,
    file_num: int,
    partitions: str = None,
):
    """
    according mzML file to map the spectrum message to parquet.
    :param parquet_path: psm_parquet_path or feature_parquet_path
    :param mzml_directory: mzml file folder
    :param output_folder: Folder where the Json file will be generated
    :param file_num: reference num
    :param partitions: The field used for splitting files, multiple fields are separated by ,
    retrun: None
    """
    if partitions:
        partitions = partitions.split(",")
    generate_features_of_spectrum(parquet_path, mzml_directory, output_folder, file_num, partitions)

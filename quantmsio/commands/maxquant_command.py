from quantmsio.core.maxquant import MaxQuant
import click
from quantmsio.core.project import create_uuid_filename

@click.command(
    "convert-maxquant-psm",
    short_help="Convert psm from mzTab to parquet file in quantms io",
)
@click.option(
    "--msms_file",
    help="the msms.txt file, this will be used to extract the peptide information",
    required=True,
)
@click.option(
    "--sdrf_file",
    help="the SDRF file needed to extract some of the metadata",
    required=True,
)
@click.option(
    "--output_folder",
    help="Folder where the parquet file will be generated",
    required=True,
)
@click.option(
    "--chunksize",
    help="Read batch size",
    default=1000000,
)
@click.option(
    "--output_prefix_file",
    help="Prefix of the parquet file needed to generate the file name",
    required=False,
)
def convert_maxquant_psm(
    msms_file: str,
    sdrf_file: str,
    output_folder: str,
    chunksize: int,
    output_prefix_file: str,
):
    """
    convert mztab psm section to a parquet file. The parquet file will contain the features and the metadata.
    :param msms_file: the msms.txt file, this will be used to extract the peptide information
    :param sdrf_file: the SDRF file needed to extract some of the metadata
    :param output_folder: Folder where the Json file will be generated
    :param chunksize: Read batch size
    :param output_prefix_file: Prefix of the Json file needed to generate the file name
    """

    if msms_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    if not output_prefix_file:
        output_prefix_file = ""

    MQ = MaxQuant(sdrf_file, msms_file)
    output_path = output_folder + "/" + create_uuid_filename(output_prefix_file, ".psm.parquet")
    MQ.convert_to_parquet(output_path=output_path, chunksize=chunksize)

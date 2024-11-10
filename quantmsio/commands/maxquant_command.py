from quantmsio.core.maxquant import MaxQuant
import click
from quantmsio.core.project import create_uuid_filename


@click.command(
    "convert-maxquant-psm",
    short_help="Convert psm from msms.txt to parquet file in quantms.io",
)
@click.option(
    "--msms_file",
    help="the msms.txt file, this will be used to extract the peptide information",
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
    output_folder: str,
    chunksize: int,
    output_prefix_file: str,
):
    """
    convert maxquant psm section to a parquet file.
    :param msms_file: the msms.txt file, this will be used to extract the peptide information
    :param output_folder: Folder where the Json file will be generated
    :param chunksize: Read batch size
    :param output_prefix_file: Prefix of the Json file needed to generate the file name
    """

    if msms_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    if not output_prefix_file:
        output_prefix_file = ""

    MQ = MaxQuant()
    output_path = output_folder + "/" + create_uuid_filename(output_prefix_file, ".psm.parquet")
    MQ.convert_psm_to_parquet(msms_path=msms_file, output_path=output_path, chunksize=chunksize)


@click.command(
    "convert-maxquant-feature",
    short_help="Convert feature from evidence to parquet file in quantms.io",
)
@click.option(
    "--evidence_file",
    help="the evidence.txt file, this will be used to extract the peptide information",
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
def convert_maxquant_feature(
    evidence_file: str,
    sdrf_file: str,
    output_folder: str,
    chunksize: int,
    output_prefix_file: str,
):
    """
    convert mztab psm section to a parquet file. The parquet file will contain the features and the metadata.
    :param evidence_file: the msms.txt file, this will be used to extract the peptide information
    :param sdrf_file: the SDRF file needed to extract some of the metadata
    :param output_folder: Folder where the Json file will be generated
    :param chunksize: Read batch size
    :param output_prefix_file: Prefix of the Json file needed to generate the file name
    """

    if evidence_file is None or sdrf_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    if not output_prefix_file:
        output_prefix_file = ""

    MQ = MaxQuant()
    output_path = output_folder + "/" + create_uuid_filename(output_prefix_file, ".psm.parquet")
    MQ.convert_feature_to_parquet(
        evidence_path=evidence_file, sdrf_path=sdrf_file, output_path=output_path, chunksize=chunksize
    )

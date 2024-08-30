import click
from quantmsio.core.maxquant_convert import MaxquantConvert
from quantmsio.core.project import create_uuid_filename

def extract_zip_list(file):
    f = open(file, encoding="utf-8")
    zip_list = f.read().splitlines()
    return zip_list

@click.command("convert-maxquant", short_help="Convert msstats/mztab to parquet file")
@click.option(
    "--sdrf_file",
    help="the SDRF file needed to extract some of the metadata",
    required=True,
)
@click.option(
    "--evidence_file",
    help="the evidence file from maxquant",
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
def convert_maxquant(
    sdrf_file: str,
    evidence_file: str,
    output_folder: str,
    output_prefix_file: str = None,
):
    """
    convert maxquant evidence to a parquet file. The parquet file will contain the features and the metadata.
    :param sdrf_file: the SDRF file needed to extract some of the metadata
    :param evidence_file: the evidence file from maxquant
    :param output_folder: Folder where the Json file will be generated
    :param output_prefix_file: Prefix of the Json file needed to generate the file name
    :return: none
    """

    if evidence_file is None or output_folder is None or sdrf_file is None:
        raise click.UsageError("Please provide all the required parameters")

    if not output_prefix_file:
        output_prefix_file = ""

    maxquant_manager = MaxquantConvert()
    output_path = output_folder + "/" + create_uuid_filename(output_prefix_file, ".maxquant.parquet")
    maxquant_manager.convert_to_parquet(evidence_file,sdrf_file,output_path,chunksize=1000000)


@click.command("convert-zip-maxquant", short_help="Convert msstats/mztab to parquet file")
@click.option(
    "--sdrf_file",
    help="the SDRF file needed to extract some of the metadata",
    required=True,
)
@click.option(
    "--maxquant_zip_files",
    help="Zip file that meets specific requirements",
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
def convert_zip_maxquant(
    sdrf_file: str,
    maxquant_zip_files: str,
    output_folder: str,
    output_prefix_file: str = None,
):
    """
    convert maxquant evidence to a parquet file. The parquet file will contain the features and the metadata.
    :param sdrf_file: the SDRF file needed to extract some of the metadata
    :param maxquant_zip_files: Zip file that meets specific requirements
    :param output_folder: Folder where the Json file will be generated
    :param output_prefix_file: Prefix of the Json file needed to generate the file name
    :return: none
    """

    if maxquant_zip_files is None or output_folder is None or sdrf_file is None:
        raise click.UsageError("Please provide all the required parameters")

    if not output_prefix_file:
        output_prefix_file = ""
    
    maxquant_manager = MaxquantConvert()
    output_path = output_folder + "/" + create_uuid_filename(output_prefix_file, ".maxquant.parquet")
    zip_list = extract_zip_list(maxquant_zip_files)
    maxquant_manager.convert_zip_to_parquet(zip_list,sdrf_file,output_path)
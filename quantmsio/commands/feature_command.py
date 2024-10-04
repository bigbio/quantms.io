import click

from quantmsio.core.feature import Feature
from quantmsio.core.project import create_uuid_filename

@click.command("convert-feature", short_help="Convert msstats/mztab to parquet file")
@click.option(
    "--sdrf_file",
    help="the SDRF file needed to extract some of the metadata",
    required=True,
)
@click.option(
    "--msstats_file",
    help="the MSstats input file, this will be considered the main format to convert",
    required=True,
)
@click.option(
    "--mztab_file",
    help="the mzTab file, this will be used to extract the protein information",
    required=True,
)
@click.option(
    "--chunksize",
    help="Read batch size",
    default=1000000,
)
@click.option(
    "--protein_file",
    help="Protein file that meets specific requirements",
    required=False,
)
@click.option(
    "--output_folder",
    help="Folder where the Json file will be generated",
    required=True,
)
@click.option(
    "--output_prefix_file",
    help="Prefix of the Json file needed to generate the file name",
    required=False,
)
def convert_feature_file(
    sdrf_file: str,
    msstats_file: str,
    mztab_file: str,
    chunksize: int,
    protein_file: str,
    output_folder: str,
    output_prefix_file: str,
):
    """
    Convert a msstats/mztab file to a parquet file. The parquet file will contain the features and the metadata.
    :param sdrf_file: the SDRF file needed to extract some of the metadata
    :param msstats_file: the MSstats input file, this will be considered the main format to convert
    :param mztab_file: the mzTab file, this will be used to extract the protein
    :param chunksize: Read batch size
    :param output_folder: Folder where the Json file will be generated
    :param output_prefix_file: Prefix of the Json file needed to generate the file name
    """

    if sdrf_file is None or msstats_file is None or mztab_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    feature_manager = Feature(mzTab_path=mztab_file,sdrf_path=sdrf_file,msstats_in_path=msstats_file)
    if not output_prefix_file:
        output_prefix_file = ""
    output_path = output_folder + "/" + create_uuid_filename(output_prefix_file, ".feature.parquet")
 
    feature_manager.write_feature_to_file(
        output_path=output_path,
        chunksize=chunksize,
        protein_file=protein_file
    )


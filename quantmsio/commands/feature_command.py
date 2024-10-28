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
    "--file_num",
    help="Read batch size",
    default=50,
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
    "--partitions",
    help="The field used for splitting files, multiple fields are separated by ,",
    required=False,
)
@click.option(
    "--output_prefix_file",
    help="Prefix of the Json file needed to generate the file name",
    required=False,
)
@click.option(
    "--duckdb_max_memory", help="The maximum amount of memory allocated by the DuckDB engine (e.g 4GB)", required=False
)
@click.option("--duckdb_threads", help="The number of threads for the DuckDB engine (e.g 4)", required=False)
def convert_feature_file(
    sdrf_file: str,
    msstats_file: str,
    mztab_file: str,
    file_num: int,
    protein_file: str,
    output_folder: str,
    partitions: str,
    output_prefix_file: str,
    duckdb_max_memory: str,
    duckdb_threads: int,
):
    """
    Convert a msstats/mztab file to a parquet file. The parquet file will contain the features and the metadata.
    :param sdrf_file: the SDRF file needed to extract some of the metadata
    :param msstats_file: the MSstats input file, this will be considered the main format to convert
    :param mztab_file: the mzTab file, this will be used to extract the protein
    :param file_num: Read batch size
    :param output_folder: Folder where the Json file will be generated
    :param partitions: The field used for splitting files, multiple fields are separated by ,
    :param output_prefix_file: Prefix of the Json file needed to generate the file name
    :param duckdb_max_memory: The maximum amount of memory allocated by the DuckDB engine (e.g 4GB)
    :param duckdb_threads: The number of threads for the DuckDB engine (e.g 4)
    """

    if sdrf_file is None or msstats_file is None or mztab_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")
    feature_manager = Feature(mzTab_path=mztab_file, sdrf_path=sdrf_file, msstats_in_path=msstats_file)
    if not output_prefix_file:
        output_prefix_file = ""
    filename = create_uuid_filename(output_prefix_file, ".feature.parquet")
    output_path = output_folder + "/" + filename
    if not partitions:
        feature_manager.write_feature_to_file(
            output_path=output_path,
            file_num=file_num,
            protein_file=protein_file,
            duckdb_max_memory=duckdb_max_memory,
            duckdb_threads=duckdb_threads,
        )
    else:
        partitions = partitions.split(",")
        feature_manager.write_features_to_file(
            output_folder=output_folder,
            filename=filename,
            partitions=partitions,
            file_num=file_num,
            protein_file=protein_file,
            duckdb_max_memory=duckdb_max_memory,
            duckdb_threads=duckdb_threads,
        )

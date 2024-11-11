import os
import click
from quantmsio.core.diann import DiaNNConvert
from quantmsio.core.project import create_uuid_filename


@click.command(
    "convert-diann",
    short_help="Convert diann_report to parquet and psm file of quantms.io format",
)
@click.option(
    "--report_path",
    help="the diann report file path",
    required=True,
)
@click.option("--qvalue_threshold", help="qvalue_threshold", required=True, default=0.05)
@click.option(
    "--mzml_info_folder",
    help="the foldef of mzml_info tsv file",
    required=True,
)
@click.option(
    "--sdrf_path",
    help="the SDRF file needed to extract some of the metadata",
    required=True,
)
@click.option(
    "--output_folder",
    help="Folder where the Json file will be generated",
    required=True,
)
@click.option(
    "--protein_file",
    help="Protein file that meets specific requirements",
    required=False,
)
@click.option(
    "--output_prefix_file",
    help="Prefix of the Json file needed to generate the file name",
    required=False,
)
@click.option(
    "--partitions",
    help="The field used for splitting files, multiple fields are separated by ,",
    required=False,
)
@click.option(
    "--duckdb_max_memory", help="The maximum amount of memory allocated by the DuckDB engine (e.g 4GB)", required=False
)
@click.option("--duckdb_threads", help="The number of threads for the DuckDB engine (e.g 4)", required=False)
@click.option(
    "--file_num",
    help="The number of files being processed at the same time",
    default=100,
)
def diann_convert_to_parquet(
    report_path: str,
    qvalue_threshold: float,
    mzml_info_folder: str,
    sdrf_path: str,
    output_folder: str,
    protein_file: str,
    output_prefix_file: str,
    partitions: str,
    duckdb_max_memory: str,
    duckdb_threads: int,
    file_num: int,
):
    """
    report_path: diann report file path
    qvalue_threshold: qvalue threshold
    mzml_info_folder: mzml info file folder
    sdrf_path: sdrf file path
    output_folder: Folder where the Json file will be generated
    param partitions: The field used for splitting files, multiple fields are separated by ,
    output_prefix_file: Prefix of the Json file needed to generate the file name
    duckdb_max_memory: The maximum amount of memory allocated by the DuckDB engine (e.g 4GB)
    duckdb_threads: The number of threads for the DuckDB engine (e.g 4)
    file_num: The number of files being processed at the same time
    """
    if report_path is None or mzml_info_folder is None or output_folder is None or sdrf_path is None:
        raise click.UsageError("Please provide all the required parameters")

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if not output_prefix_file:
        output_prefix_file = ""
    filename = create_uuid_filename(output_prefix_file, ".feature.parquet")
    feature_output_path = output_folder + "/" + filename

    dia_nn = DiaNNConvert(
        diann_report=report_path,
        sdrf_path=sdrf_path,
        duckdb_max_memory=duckdb_max_memory,
        duckdb_threads=duckdb_threads,
    )
    if not partitions:
        dia_nn.write_feature_to_file(
            qvalue_threshold=qvalue_threshold,
            mzml_info_folder=mzml_info_folder,
            output_path=feature_output_path,
            file_num=file_num,
            protein_file=protein_file,
        )
    else:
        partitions = partitions.split(",")
        dia_nn.write_features_to_file(
            qvalue_threshold=qvalue_threshold,
            mzml_info_folder=mzml_info_folder,
            output_folder=output_folder,
            filename=filename,
            partitions=partitions,
            file_num=file_num,
            protein_file=protein_file,
        )

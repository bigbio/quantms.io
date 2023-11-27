from quantms_io.core.diann_convert import DiaNNConvert
import click
from quantms_io.core.project import create_uuid_filename
from typing import List
import os

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    This is the main tool that gives access to all commands.
    """

@click.command("convert-diann", short_help="Convert diann_report to parquet and psm file of quantms.io format")
@click.option(
    "--report_path",
    help="the diann report file path",
    required=True,
)
@click.option(
    "--design_file",
    help="the disign file path",
    required=True,
)
@click.option(
    "--modifications",
    nargs=2, 
    type=str,
    help="a list contains fix modifications and variable modifications",
    required=True,
)
@click.option(
    "--qvalue_threshold",
    help="qvalue_threshold",
    required=True,
    default= 0.05
)
@click.option(
    "--mzml_info_folder",
    help="mzml info tsv file",
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
@click.option("--output_prefix_file", help="Prefix of the Json file needed to generate the file name", required=False)
@click.option("--duckdb_max_memory", help= "The maximum amount of memory allocated by the DuckDB engine (e.g 4GB)")
@click.option("--duckdb_threads", help= "The number of threads for the DuckDB engine (e.g 4)")
@click.option("--file_num", help= "The number of files being processed at the same time", default = 100)
@click.pass_context
def diann_convert_to_parquet(ctx, report_path: str, design_file: str, modifications:List, qvalue_threshold: float,
                             mzml_info_folder:str, sdrf_path:str, output_folder:str, output_prefix_file:str,
                             duckdb_max_memory:str, duckdb_threads:int, file_num:int ):
    '''
    report_path: diann report file path
    design_file: the disign file path
    modifications: a list contains fix modifications and variable modifications
    qvalue_threshold: qvalue threshold
    mzml_info_folder: mzml info file folder
    sdrf_path: sdrf file path
    output_folder: Folder where the Json file will be generated
    output_prefix_file: Prefix of the Json file needed to generate the file name
    duckdb_max_memory: The maximum amount of memory allocated by the DuckDB engine (e.g 4GB)
    duckdb_threads: The number of threads for the DuckDB engine (e.g 4)
    file_num: The number of files being processed at the same time
    '''
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    if not output_prefix_file:
        output_prefix_file = ''
    
    feature_output_path = output_folder + "/" + create_uuid_filename(output_prefix_file,'.feature.parquet')
    psm_output_path = output_folder + "/" + create_uuid_filename(output_prefix_file,'.psm.parquet')

    DiaNN = DiaNNConvert()

    DiaNN.generate_psm_and_feature_file(
                                        report_path=report_path,
                                        qvalue_threshold=qvalue_threshold,
                                        mzml_info_folder=mzml_info_folder,
                                        design_file=design_file,
                                        modifications=modifications,
                                        sdrf_path = sdrf_path,
                                        psm_output_path=psm_output_path,
                                        feature_output_path = feature_output_path,
                                        duckdb_max_memory= duckdb_max_memory,
                                        duckdb_threads= duckdb_threads,
                                        file_num = file_num
                                    )

cli.add_command(diann_convert_to_parquet)
if __name__ == '__main__':
    cli()
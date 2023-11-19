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

@click.command(
    "diann_convert_to_parquet", short_help="Convert diann_report to parquet and psm file"
)
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
    "--fasta_path",
    help="reference fasta database file",
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
    "--pg_path",
    help="pg_matrix file",
    required=True,
)
@click.option(
    "--pr_path",
    help="pr_matrix file",
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
    help="sdrf file path",
    required=True,
)
@click.option(
    "--output_folder",
    help="Folder where the Json file will be generated",
    required=True,
)
@click.option("--output_prefix_file", help="Prefix of the Json file needed to generate the file name", required=False)
@click.option(
    "--chunksize",
    help="mzml info tsv file",
    required=False,
    default = 100000
)
@click.pass_context
def diann_convert_to_parquet(ctx,report_path:str,design_file:str,fasta_path:str,modifications:List,pg_path:str,pr_path:str,qvalue_threshold: float,mzml_info_folder:str,sdrf_path:str,output_folder:str,output_prefix_file:str,chunksize:int):
    '''
    report_path: diann report file path
    design_file: the disign file path
    fasta_path: reference fasta database file
    modifications: a list contains fix modifications and variable modifications
    pg_path: pg_matrix file
    pr_path: pr_matrix file
    qvalue_threshold: qvalue threshold
    mzml_info_folder: mzml info file folder
    sdrf_path: sdrf file path
    output_folder: Folder where the Json file will be generated
    output_prefix_file: Prefix of the Json file needed to generate the file name
    chunksize: batch size
    '''
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    if not output_prefix_file:
        output_prefix_file = ''
    
    feature_output_path = output_folder + "/" + create_uuid_filename(output_prefix_file,'.feature.parquet')
    psm_output_path = output_folder + "/" + create_uuid_filename(output_prefix_file,'.psm.parquet')

    DiaNN = DiaNNConvert()

    DiaNN.generate_feature_and_psm_file(report_path=report_path,design_file=design_file,fasta_path=fasta_path,
                                        modifications=modifications,pg_path=pg_path,pr_path=pr_path,
                                        qvalue_threshold=qvalue_threshold,mzml_info_folder=mzml_info_folder,
                                        sdrf_path=sdrf_path,output_path=feature_output_path,
                                        psm_output_path=psm_output_path,chunksize=chunksize)
    #DiaNN.generate_psm_file(report_path,design_file,fasta_path,modifications,pg_path,qvalue_threshold,mzml_info_folder,psm_output_path,chunksize)

cli.add_command(diann_convert_to_parquet)
if __name__ == '__main__':
    cli()
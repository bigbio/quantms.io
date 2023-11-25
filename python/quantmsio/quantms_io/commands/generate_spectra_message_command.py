from quantms_io.core.tools import generate_features_of_spectrum
import click
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    This is the main tool that gives access to all commands.
    """

@click.command("map_spectrum_message_to_parquet", short_help="According mzMl to map the spectrum message to parquet"
)
@click.option('--parquet_path',  help='Psm or feature parquet path')
@click.option('--mzml_directory', help='mzml file folder')
@click.option('--output_path', help='output file path(.parquet)')
@click.option('--label', type=click.Choice(['feature', 'psm'], case_sensitive=False),help='parquet type')
@click.option("--chunksize",help="The number of rows of parquet read using pandas streaming", default=100000)
@click.option('--partition', type=click.Choice(['charge', 'reference_file_name'], case_sensitive=False),help='partition',required=False)
@click.pass_context
def map_spectrum_message_to_parquet(ctx,parquet_path: str, mzml_directory: str,output_path:str,label:str,chunksize,partition:str=None):
    '''
    according mzML file to map the spectrum message to parquet.
    :param parquet_path: psm_parquet_path or feature_parquet_path
    :param mzml_directory: mzml file folder
    :param output_path: output file path(extension[.parquet])
    :param label: feature or psm
    :param chunksize: batch size
    :param partition: charge or reference_file_name
    retrun: None
    '''
    if not output_path.endswith('parquet'):
        raise click.UsageError("Please provide file extension(.parquet)")
    generate_features_of_spectrum(
        parquet_path,
        mzml_directory,
        output_path,
        label,
        chunksize,
        partition
    )

cli.add_command(map_spectrum_message_to_parquet)
if __name__ == '__main__':
    cli()
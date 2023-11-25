import click

from quantms_io.core.psm import PSMHandler
from quantms_io.core.project import create_uuid_filename
from quantms_io.core.tools import plot_peptidoform_charge_venn, plot_sequence_venn
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    This is the main tool that gives access to all commands.
    """

@click.command("convert-psm", short_help="Convert psm from mzTab to parquet file in quantms io",)
@click.option(
    "--mztab_file",
    help="the mzTab file, this will be used to extract the protein information",
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
              required=False)
@click.option(
    "--verbose",
              help="Output debug information.",
    default=False, is_flag=True)
@click.pass_context
def convert_psm_file(ctx, mztab_file: str, output_folder: str, output_prefix_file: str=None, verbose: bool = False):
    """
    convert mztab psm section to a parquet file. The parquet file will contain the features and the metadata.
    :param mztab_file: the mzTab file, this will be used to extract the protein information
    :param output_folder: Folder where the Json file will be generated
    :param output_prefix_file: Prefix of the Json file needed to generate the file name
    :param verbose: Output debug information.
    :return: none
    """

    if mztab_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")


    if not output_prefix_file:
        output_prefix_file = ''

    psm_manager = PSMHandler()
    psm_manager.parquet_path = output_folder + "/" + create_uuid_filename(output_prefix_file,'.psm.parquet')
    psm_manager.convert_mztab_to_psm(
        mztab_path=mztab_file, output_folder=output_folder, parquet_path=psm_manager.parquet_path, verbose=verbose
    )


@click.command(
    "compare_set_of_psms", short_help="plot venn for a set of Psms parquet"
)
@click.option('-p','--parquets', type=str, help='List of psm parquet path', multiple=True)
@click.option('-t','--tags', type=str, help='List of parquet label', multiple=True)
@click.pass_context
def compare_set_of_psms(ctx,parquets, tags):
    """
    Compare a set of psm parquet files
    :param parquets: a set of psm parquet path
    :param tags: a set of psm label
    """
    if len(parquets) != len(tags):
        raise click.UsageError("Please provide same length of parquet_list and label_list")

    plot_peptidoform_charge_venn(parquets, tags)
    plot_sequence_venn(parquets, tags)

cli.add_command(convert_psm_file)
cli.add_command(compare_set_of_psms)

if __name__ == '__main__':
    cli()
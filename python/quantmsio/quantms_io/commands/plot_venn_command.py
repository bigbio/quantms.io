import click
from quantms_io.core.tools import plot_peptidoform_charge_venn,plot_sequence_venn
@click.command(
    "get message from a set of PSMS", short_help="plot venn for a set of Psms"
)
@click.option('--parquets', nargs=-1, type=str, help='List of psm parquet path')
@click.option('--tags', nargs=-1, type=str, help='List of parquet label')
def compare_set_of_psms(parquets,tags):
    '''
    :param parquet_list: a set of psm parquet path
    :param label_list: a set of psm label
    '''
    if len(parquets) != len(tags):
        raise click.UsageError("Please provide same length of parquet_list and label_list")
    
    plot_peptidoform_charge_venn(parquets,tags)
    plot_sequence_venn(parquets,tags)

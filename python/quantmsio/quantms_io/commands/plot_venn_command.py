import click
from quantms_io.core.tools import plot_peptidoform_charge_venn,plot_sequence_venn
@click.command(
    "get message from a set of PSMS", short_help="Convert msstats/mztab to parquet file"
)
@click.option('--parquet_list', nargs=-1, type=str, help='List of psm parquet path')
@click.option('--label_list', nargs=-1, type=str, help='List of parquet label')
def compare_set_of_psms(parquet_list,label_list):
    '''
    :param parquet_list: a set of psm parquet path
    :param label_list: a set of psm label
    '''
    if len(parquet_list) != len(label_list):
        raise click.UsageError("Please provide same length of parquet_list and label_list")
    
    plot_peptidoform_charge_venn(parquet_list,label_list)
    plot_sequence_venn(parquet_list,label_list)

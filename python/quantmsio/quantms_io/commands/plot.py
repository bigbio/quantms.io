import click
from quantms_io.core.tools import plot_peptides_of_lfq_condition
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    This is the main tool that gives access to all commands.
    """

@click.command("plot_peptides", short_help="plot peptides of condition in lfq",)
@click.option(
    "--psm_parquet_path",
    help="psm parquet path in lfq",
    required=True,
)
@click.option(
    "--sdrf_path",
    help="sdrf path",
    required=True,
)
@click.option(
    "--save_path",
    help="img save path [xxx.png]",
    required=True)
@click.pass_context
def plot_peptides(ctx, psm_parquet_path: str, sdrf_path: str, save_path: str):
    """
    convert mztab psm section to a parquet file. The parquet file will contain the features and the metadata.
    :param psm_parquet_path: psm parquet path in lfq
    :param sdrf_path: sdrf path
    :param save_path: img save path [xxx.png]
    :return: none
    """
    plot_peptides_of_lfq_condition(psm_parquet_path=psm_parquet_path,sdrf_path=sdrf_path,save_path=save_path)


cli.add_command(plot_peptides)

if __name__ == '__main__':
    cli()
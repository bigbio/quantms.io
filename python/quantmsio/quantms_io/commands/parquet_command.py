import click
import pandas as pd
import datacompy
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    This is the main tool that gives access to all commands.
    """

@click.command("compare-parquet", short_help="Compare two parquet files, feature or psm files")
@click.option("--parquet_path_one", help="First parquet file", required=True)
@click.option("--parquet_path_two", help="the parquet file of memory version", required=True)
@click.option("--report_path", help="report path", required=True)
@click.pass_context
def compare_two_parquet(ctx,parquet_path_one: str, parquet_path_two: str, report_path: str):
    """
    compare two parquet files generated from the same mztab psms files.
    :param parquet_path_one: the parquet file of discache version
    :param parquet_path_two: the parquet file of memory version
    :param report_path: report path
    """
    parquet_one = pd.read_parquet(parquet_path_one)
    parquet_two = pd.read_parquet(parquet_path_two)
    parquet_one = parquet_one.astype(str)
    parquet_two = parquet_two.astype(str)
    compare = datacompy.Compare(parquet_one, parquet_two, join_columns='sequence', df1_name='discache',df2_name='no_cache')
    with open(report_path,'w') as f:
        f.write(compare.report())

cli.add_command(compare_two_parquet)
if __name__ == '__main__':
    cli()
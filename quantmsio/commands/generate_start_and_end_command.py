import click

from quantmsio.core.tools import generate_start_and_end_from_fasta


@click.command(
    "inject-start-and-end-from-fasta",
    short_help="inject start and end to feature or psm",
)
@click.option("--parquet_path", help="psm or feature parquet file path", required=True)
@click.option("--fasta_path", help="diann feature parquet file path", required=True)
@click.option(
    "--label",
    type=click.Choice(["feature", "psm"], case_sensitive=False),
    help="parquet type",
)
@click.option("--output_path", help="save path", required=True)
def inject_start_and_end_from_fasta(parquet_path: str, fasta_path: str, label: str, output_path: str):
    """
    Register the file with project.json
    :param parquet_path: psm or feature parquet file path
    :param fasta_path: fasta file path
    :param label parquet file type (psm or feature)
    :param output_path: save path
    :return: none
    """
    generate_start_and_end_from_fasta(
        parquet_path=parquet_path,
        fasta_path=fasta_path,
        label=label,
        output_path=output_path,
    )

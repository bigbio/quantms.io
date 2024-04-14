import click

from quantmsio.core.tools import map_protein_for_parquet
from quantmsio.core.tools import map_protein_for_tsv

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


# Command Group
@click.group(name="labels", context_settings=CONTEXT_SETTINGS)
def labels():
    """Tool related commands"""
    pass


@labels.command(
    "convert-accession",
    short_help="Map reported proteins to uniprot names/accessions in fasta file",
)
@click.option("--parquet_path", help="psm or feature parquet path")
@click.option("--fasta", help="Reference fasta database")
@click.option("--output_path", help="output file path")
@click.option(
    "--map_parameter",
    type=click.Choice(["map_protein_name", "map_protein_accession"], case_sensitive=False),
    help="map type",
)
@click.option(
    "--label",
    type=click.Choice(["feature", "psm"], case_sensitive=False),
    help="parquet type",
)
def get_unanimous_for_parquet(parquet_path, fasta, output_path, map_parameter, label):
    """
    according fasta database to map the proteins accessions to uniprot names.
    :param parquet_path: psm_parquet_path or feature_parquet_path
    :param fasta: Reference fasta database
    :param output_path: output file path
    :param map_parameter: map_protein_name or map_protein_accession
    :param label: feature or psm
    return: None
    """
    map_protein_for_parquet(parquet_path, fasta, output_path, map_parameter, label)


# tsv
@labels.command(
    "get-unanimous-for-tsv",
    short_help="According fasta database to map the proteins accessions to uniprot names.",
)
@click.option("--path", help="ae or de path")
@click.option("--fasta", help="Reference fasta database")
@click.option("--output_path", help="output file path")
@click.option(
    "--map_parameter",
    type=click.Choice(["map_protein_name", "map_protein_accession"], case_sensitive=False),
    help="map type",
)
def get_unanimous_for_tsv(path, fasta, output_path, map_parameter):
    """
    according fasta database, to map the proteins accessions to uniprot names.
    :param path: de_path or ae_path
    :param fasta: Reference fasta database
    :param output_path: output file path
    :param map_parameter: map_protein_name or map_protein_accession
    retrun: None
    """
    map_protein_for_tsv(path, fasta, output_path, map_parameter)

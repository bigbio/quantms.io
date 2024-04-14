import click

from quantmsio.core.tools import map_gene_msgs_to_parquet


@click.command("map-gene-msg-to-parquet", short_help="map the gene message to parquet")
@click.option("--parquet_path", help="Psm or feature parquet path")
@click.option("--fasta_path", help="fasta file path")
@click.option("--output_path", help="output file path(.parquet)")
@click.option(
    "--label",
    type=click.Choice(["feature", "psm"], case_sensitive=False),
    help="parquet type",
)
@click.option(
    "--map_parameter",
    type=click.Choice(["map_protein_name", "map_protein_accession"], case_sensitive=False),
    help="map type",
)
@click.option("--species", help="species", default="human")
def map_gene_msg_to_parquet(
    parquet_path: str,
    fasta_path: str,
    output_path: str,
    label: str,
    map_parameter: str,
    species: str,
):
    """
    according mzML file to map the spectrum message to parquet.
    :param parquet_path: psm_parquet_path or feature_parquet_path
    :param fasta_path: fasta file path
    :param output_path: output file path(extension[.parquet])
    :param label: feature or psm
    :param map_parameter: map_protein_name or map_protein_accession
    retrun: None
    """
    if not output_path.endswith("parquet"):
        raise click.UsageError("Please provide file extension(.parquet)")

    map_gene_msgs_to_parquet(parquet_path, fasta_path, map_parameter, output_path, label, species)

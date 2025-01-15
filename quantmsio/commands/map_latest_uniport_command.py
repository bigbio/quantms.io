import click
from quantmsio.core.project import create_uuid_filename
from quantmsio.operate.tools import map_peptide_to_protein


@click.command(
    "map-latest-uniport",
    short_help="Map the peptides to the latest UniProt Fasta file.",
)
@click.option(
    "--feature_file",
    help="feature file",
    required=True,
)
@click.option(
    "--fasta",
    help="the latest UniProt Fasta file",
    required=True,
)
@click.option(
    "--output_folder",
    help="Folder where the Json file will be generated",
    required=True,
)
@click.option(
    "--output_prefix_file",
    help="Prefix of the parquet file needed to generate the file name",
    required=False,
)
def map_latest_uniport(
    feature_file: str,
    fasta: str,
    output_folder: str,
    output_prefix_file: str,
):
    """
    :param feature_file: feature file
    :param sdrf_file: the SDRF file needed to extract some of the metadata
    :param output_folder: Folder where the Json file will be generated
    :param output_prefix_file: Prefix of the Json file needed to generate the file name
    """

    if feature_file is None or fasta is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    if not output_prefix_file:
        output_prefix_file = "feature"

    filename = create_uuid_filename(output_prefix_file, ".feature.parquet")
    map_peptide_to_protein(feature_file, fasta, output_folder, filename)

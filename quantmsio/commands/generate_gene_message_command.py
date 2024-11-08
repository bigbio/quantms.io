import click
from quantmsio.operate.tools import generate_feature_of_gene


@click.command(
    "map-gene-message-to-parquet",
    short_help="According mzMl to map the gene message to parquet",
)
@click.option("--parquet_path", help="Psm or feature parquet path")
@click.option("--fasta", help="fasta file")
@click.option(
    "--output_folder",
    help="Folder where the Json file will be generated",
    required=True,
)
@click.option(
    "--file_num",
    help="The number of rows of parquet read using pandas streaming",
    default=10,
)
@click.option(
    "--partitions",
    help="The field used for splitting files, multiple fields are separated by ,",
    required=False,
)
def map_gene_message_to_parquet(
    parquet_path: str,
    fasta: str,
    output_folder: str,
    file_num: int,
    partitions: str = None,
):
    """
    according fasta file to map the gene message to parquet.
    :param parquet_path: psm_parquet_path or feature_parquet_path
    :param fasta: fasta path
    :param output_folder: Folder where the Json file will be generated
    :param file_num: reference num
    :param partitions: The field used for splitting files, multiple fields are separated by ,
    retrun: None
    """
    if partitions:
        partitions = partitions.split(",")
    generate_feature_of_gene(parquet_path, fasta, output_folder, file_num, partitions)
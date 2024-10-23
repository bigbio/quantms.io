import click
from quantmsio.core.project import create_uuid_filename
from quantmsio.core.psm import Psm
from quantmsio.operate.plots import plot_peptidoform_charge_venn, plot_sequence_venn


@click.command(
    "convert-psm",
    short_help="Convert psm from mzTab to parquet file in quantms io",
)
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
    "--chunksize",
    help="Read batch size",
    default=1000000,
)
@click.option(
    "--protein_file",
    help="Protein file that meets specific requirements",
    required=False,
)
@click.option(
    "--output_prefix_file",
    help="Prefix of the parquet file needed to generate the file name",
    required=False,
)
def convert_psm_file(
    mztab_file: str,
    output_folder: str,
    chunksize: int,
    protein_file: str,
    output_prefix_file: str,
):
    """
    convert mztab psm section to a parquet file. The parquet file will contain the features and the metadata.
    :param mztab_file: the mzTab file, this will be used to extract the protein information
    :param output_folder: Folder where the Json file will be generated
    :param chunksize: Read batch size
    :param output_prefix_file: Prefix of the Json file needed to generate the file name
    """

    if mztab_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    if not output_prefix_file:
        output_prefix_file = ""

    psm_manager = Psm(mzTab_path=mztab_file)
    output_path = output_folder + "/" + create_uuid_filename(output_prefix_file, ".psm.parquet")
    psm_manager.write_psm_to_file(output_path=output_path, chunksize=chunksize, protein_file=protein_file)


@click.command("compare-set-psms", short_help="plot venn for a set of Psms parquet")
@click.option("-p", "--parquets", type=str, help="List of psm parquet path", multiple=True)
@click.option("-t", "--tags", type=str, help="List of parquet label", multiple=True)
def compare_set_of_psms(parquets, tags):
    """
    Compare a set of psm parquet files
    :param parquets: a set of psm parquet path
    :param tags: a set of psm label
    """
    if len(parquets) != len(tags):
        raise click.UsageError("Please provide same length of parquet_list and label_list")

    plot_peptidoform_charge_venn(parquets, tags)
    plot_sequence_venn(parquets, tags)

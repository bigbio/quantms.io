import click

from quantmsio.core.project import create_uuid_filename
from quantmsio.core.psm import PSMHandler
from quantmsio.core.tools import plot_peptidoform_charge_venn
from quantmsio.core.tools import plot_sequence_venn


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
    "--use_cache",
    help="Use cache instead of in memory conversion",
    required=False,
    is_flag=True,
)
@click.option(
    "--output_prefix_file",
    help="Prefix of the parquet file needed to generate the file name",
    required=False,
)
@click.option("--verbose", help="Output debug information.", default=False, is_flag=True)
def convert_psm_file(
    mztab_file: str,
    output_folder: str,
    use_cache: bool,
    output_prefix_file: str = None,
    verbose: bool = False,
):
    """
    convert mztab psm section to a parquet file. The parquet file will contain the features and the metadata.
    :param mztab_file: the mzTab file, this will be used to extract the protein information
    :param output_folder: Folder where the Json file will be generated
    :param output_prefix_file: Prefix of the Json file needed to generate the file name
    :param verbose: Output debug information.
    :return: none
    """

    if mztab_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    if not output_prefix_file:
        output_prefix_file = ""

    psm_manager = PSMHandler()
    psm_manager.parquet_path = output_folder + "/" + create_uuid_filename(output_prefix_file, ".psm.parquet")
    psm_manager.convert_mztab_to_psm(
        mztab_path=mztab_file,
        parquet_path=psm_manager.parquet_path,
        verbose=verbose,
        use_cache=use_cache,
    )


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

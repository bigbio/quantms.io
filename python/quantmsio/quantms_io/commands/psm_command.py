import click

from quantms_io.core.feature import FeatureHandler
from quantms_io.core.psm import PSMHandler


@click.command(
    "convert-psm-file",
    short_help="Convert psm from mzTab to parquet file in quantms io",
)
@click.option(
    "--mztab_file",
    help="the mzTab file, this will be used to extract the protein information",
    required=True,
)
@click.option(
    "--output_folder",
    help="Folder where the Json file will be generated",
    required=True,
)
@click.option("--output_file", help="Parquet file name", required=False)
def convert_psm_file(mztab_file: str, output_folder: str, output_file: str):
    """
    Convert mztab psm section to a parquet file. The parquet file will contain the features and the metadata.
    :param mztab_file: the mzTab file, this will be used to extract the protein information
    :param output_folder: Folder where the Json file will be generated
    :param output_file: Prefix of the Json file needed to generate the file name
    :return: none
    """

    if mztab_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    psm_manager = PSMHandler()
    psm_manager.parquet_path = output_folder + "/" + output_file
    psm_manager.convert_mztab_to_feature(
        mztab_path=mztab_file, parquet_path=psm_manager.parquet_path
    )

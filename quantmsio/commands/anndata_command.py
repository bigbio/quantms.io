import click
import pandas as pd
from quantmsio.core.combiner import Combiner
from quantmsio.utils.file_utils import find_ae_files
from quantmsio.core.project import create_uuid_filename


@click.command(
    "merge-ae-files",
    short_help="Merge multiple AE files into a file in AnnData format.",
)
@click.option(
    "--directory",
    help="The directory for storing AE files.",
    required=True,
)
@click.option("--output_folder", help="Folder to generate the adata file.", required=True)
@click.option("--output_prefix_file", help="Prefix of the df expression file", required=False)
def merge_ae_files(
    directory: str,
    output_folder: str,
    output_prefix_file: str,
):
    ae_files = find_ae_files(directory)
    output_path = output_folder + "/" + create_uuid_filename(output_prefix_file, ".h5ad")
    ae_combiner = Combiner()
    if len(ae_files) == 0:
        raise click.UsageError("No AE files were found.")
    else:
        for file in ae_files:
            df = pd.read_csv(file, comment="#", sep="\t")
            adata = ae_combiner.transform_to_adata(df)
            ae_combiner.combine_adata(adata)
        ae_combiner.save_adata(output_path)

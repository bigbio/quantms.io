from pathlib import Path
from typing import Optional

import click

from quantmsio.core.fragpipe import FragPipe


@click.command(
    "convert-fragpipe-psm",
    short_help="Convert FragPipe PSMs from psm.tsv to parquet file in quantms.io",
)
@click.option(
    "--msms_file",
    type=click.Path(dir_okay=False, path_type=Path, file_okay=True, exists=True),
    help="the psm.tsv file, this will be used to extract the peptide information",
    required=True,
)
@click.option(
    "-o",
    "--output_folder",
    type=click.Path(dir_okay=True, path_type=Path, file_okay=False),
    help="Folder where the parquet file will be generated",
    required=True,
)
@click.option("-b", "--chunksize", help="Read batch size", default=1000000, type=int)
@click.option(
    "--output_prefix_file",
    help="Prefix of the parquet file needed to generate the file name",
    required=False,
)
def convert_fragpipe_psm(
    msms_file: Path,
    output_folder: Path,
    chunksize: int,
    output_prefix_file: Optional[str] = None,
):
    if not output_folder.exists():
        output_folder.mkdir(parents=True, exist_ok=True)
    converter = FragPipe(output_directory=output_folder)
    converter.write_psms_to_parquet(msms_file, batch_size=chunksize, output_prefix_file=output_prefix_file)

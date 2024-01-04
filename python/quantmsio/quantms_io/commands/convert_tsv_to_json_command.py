from quantms_io.core.json import JsonConverter
import click
import os

@click.command("convert-tsv-to-json", short_help="Convert AE or DE file to JSON format", )
@click.option("--file", help="AE or DE file", required=True)
def convert_tsv_to_json(file: str):
    if not os.path.exists(file):
        raise click.UsageError("The file does not exist.")

    converter = JsonConverter()
    converter.convert_tsv_to_json(file)



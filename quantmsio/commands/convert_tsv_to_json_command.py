import os

import click

from quantmsio.core.json import JsonConverter

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


# Command Group
@click.group(name="json", context_settings=CONTEXT_SETTINGS)
def json():
    """Tool related commands"""
    pass


@json.command(
    "convert-tsv-to-json",
    short_help="Convert AE or DE file to JSON format",
)
@click.option("--file", help="AE or DE file", required=True)
@click.option("--json_path", help="save json path", required=True)
def convert_tsv_to_json(file: str, json_path: str):
    if not os.path.exists(file):
        raise click.UsageError("The file does not exist.")

    converter = JsonConverter()
    converter.convert_tsv_to_json(file, json_path)


@json.command(
    "convert-sdrf-to-json",
    short_help="Convert sdrf file to JSON format",
)
@click.option("--file", help="sdrf file", required=True)
@click.option("--json_path", help="save json path", required=True)
def convert_sdrf_to_json(file: str, json_path: str):
    if not os.path.exists(file):
        raise click.UsageError("The file does not exist.")

    converter = JsonConverter()
    converter.sdrf_to_json(file, json_path)

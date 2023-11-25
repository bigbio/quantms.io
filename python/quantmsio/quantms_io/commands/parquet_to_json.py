import click

from quantms_io.core.json import JsonConverter

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    This is the main tool that gives access to all commands.
    """

@click.command("convert-parquet-json", short_help="Convert parquet to json")
@click.option("--data_type", type=click.Choice(['feature', 'psm']),
              help="Data type of the parquet: psm or feature", required=True)
@click.option("--parquet_path", help="The parquet path", required=True)
@click.option("--json_path", help="The json path", required=True)
@click.pass_context
def convert_parquet_to_json(ctx, data_type: str, parquet_path: str, json_path: str):
    """
    Convert parquet to json
    :param ctx:
    :param data_type:
    :param parquet_path:
    :param json_path:
    :return:
    """
    converter = JsonConverter()

    if data_type == 'psm':
        converter.psm_to_json(parquet_path, json_path)
    elif data_type == 'feature':
        converter.feature_to_json(parquet_path, json_path)
    else:
        raise ValueError(f"Unknown data type {data_type}")
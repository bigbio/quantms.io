from quantms_io.core.tools import convert_to_json
import click
import os
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    This is the main tool that gives access to all commands.
    """


@click.command("convert_tsv_to_json", short_help="Register the file to project.json.",)
@click.option("--file", help="AE or DE file", required=True)
@click.pass_context
def convert_tsv_to_json(ctx,file):
    if not os.path.exists(file):
        raise click.UsageError("The file does not exist.")
    
    convert_to_json(file)

cli.add_command(convert_tsv_to_json)
if __name__ == '__main__':
    cli()
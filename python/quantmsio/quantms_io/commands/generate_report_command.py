import click

from quantms_io.core.tools import generate_report_of_psms_or_features
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    This is the main tool that gives access to all commands.
    """

@click.command("generate-report", short_help="generate report of psm or feature files " "format",)
@click.option(
    "--check_dir", help="Folder to generate the df expression file.", required=True
)
@click.option('--label', type=click.Choice(['feature', 'psm'], case_sensitive=False),help='parquet type')
@click.pass_context
def generate_report_about_files(ctx, check_dir: str, label: str):
    '''
    ckeck_dir: psm or feature file directory
    label: psm or feature
    '''
    generate_report_of_psms_or_features(check_dir=check_dir,label=label)

cli.add_command(generate_report_about_files)

if __name__ == '__main__':
    cli()
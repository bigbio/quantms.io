import click

from quantms_io.core.tools import generate_report_of_psms_or_features

@click.command("generate-report", short_help="generate report of psm or feature files " "format",)
@click.option(
    "--check_dir", help="Folder to generate the df expression file.", required=True
)
@click.option('--label', type=click.Choice(['feature', 'psm'], case_sensitive=False),help='parquet type')
def generate_report_about_files(check_dir: str, label: str):
    """
    ckeck_dir: psm or feature file directory
    label: psm or feature
    """
    generate_report_of_psms_or_features(check_dir=check_dir,label=label)

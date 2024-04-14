import click

from quantmsio.core.tools import generate_project_report


@click.command(
    "generate-project-report",
    short_help="generate report of project " "format",
)
@click.option("--project_folder", help="Folder to generate the df expression file.", required=True)
def generate_report_about_project(project_folder):
    """
    project_folder: The folder path for the full project.
    """
    generate_project_report(project_folder)

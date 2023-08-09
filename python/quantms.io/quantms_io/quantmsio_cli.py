"""
Commandline interface for quantmsio package allows generating the quantms.io file format from different sources and
 steps. The quantms.io specification is available in the docs folder of this repository.
"""

import click

from quantms_io import __version__ as __version__

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.version_option(version=__version__, package_name="quantms_io", message="%(package)s %(version)s")
@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    This is the main tool that gives access to all commands to convert SDRF files into pipelines specific configuration files
    """
    pass


@click.command("generate-project-json", short_help="Generate a json project file from original PX accession")
@click.option("--project_accession", "-p", help="Project accession to generate the json file", required=True)
@click.option("--sdrf_file", "-s", help="the SDRF file needed to extract some of the metadata", required=True)
@click.option("--quantms_version", "-q", help="Quantms.io version to use as output in the project file", required=True)
@click.option("--output_folder", "-o", help="Folder where the Json file will be generated", required=True)
@click.option("--output_file", "-o", help="Prefix of the Json file needed to generate the file name", required=True)
def generate_project_json(project_accession: str,
                          sdrf_file: str, quantms_version: str,
                          output_folder: str,
                          output_file: str):
    """
    Generate a json project file from the original PX accession and SDRF file. The Json file definition is available in the docs
    folder of this repository https://github.com/bigbio/quantms.io/blob/main/docs/PROJECT.md. This command will generate
    the Json and attach to the project files the SDRF file provided as input.
    :param project_accession: Project accession to generate the json file
    :param sdrf_file: the SDRF file needed to extract some of the metadata
    :param quantms_version: Quantms.io version to use as output in the project file
    :param output_folder: Folder where the Json file will be generated
    :param output_file: Prefix of the Json file needed to generate the file name
    :return: none
    """

    pass


cli.add_command(generate_project_json)


def quantms_io_main():
    """
    Main function to run the quantmsio command line interface
    :return: none
    """
    cli()


if __name__ == '__main__':
    quantms_io_main()

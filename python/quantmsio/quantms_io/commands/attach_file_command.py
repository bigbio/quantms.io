from quantms_io.core.tools import register_file_to_json
import click
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    This is the main tool that gives access to all commands.
    """


@click.command("attach-file", short_help="Register the file to project.json.",)
@click.option("--project_file", help="the project.json file", required=True)
@click.option("--attach_file", help="The path of the file that will be registered", required=True)
@click.option("--category", type=click.Choice(['feature_file', 'psm_file', 'differential_file',
                                              'absolute_file'], case_sensitive=False),
              help="The type of file that will be registered.", required=True)
@click.option("--replace_existing", help="Whether to delete old files", is_flag=True)
@click.pass_context
def attach_file_to_json(ctx, project_file, attach_file, category, replace_existing):
    """
    Register the file with project.json
    :param project_file: the project.json file path
    :param attach_file: The path of the file that will be registered
    :param category: The type of file that will be registered
    :param replace_existing Whether to delete old files
    :return: none
    """
    register_file_to_json(project_file, attach_file,
                          category, replace_existing)


cli.add_command(attach_file_to_json)

if __name__ == '__main__':
    cli()

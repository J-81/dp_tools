# Connect all sub cli to a top level namespace

import click

from dp_tools.scripts.vv_interface import validation
from dp_tools.scripts.osd_api_cli import osd

@click.group()
def cli():
    pass

cli.add_command(validation)
cli.add_command(osd)

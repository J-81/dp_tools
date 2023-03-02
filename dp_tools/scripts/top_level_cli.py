# Connect all sub cli to a top level namespace

import click

from dp_tools.scripts.vv_interface import validation
from dp_tools.scripts.osd_api_cli import osd
from dp_tools.scripts.data_assets_cli import data_assets

@click.group()
def cli():
    pass

cli.add_command(validation)
cli.add_command(osd)
cli.add_command(data_assets)

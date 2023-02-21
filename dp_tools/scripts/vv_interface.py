from pathlib import Path

import click

from dp_tools.plugin_api import load_plugin
from dp_tools.core.loaders import load_data

@click.group()
def cli():
    pass


@click.group()
def validation():
    pass

cli.add_command(validation)

@click.command()
@click.option('--output', default="VV_report.tsv", help="Name of report output file", show_default=True)
@click.argument('plug_in_dir')
@click.argument('data_dir')
def run(plug_in_dir, output):
    plugin = load_plugin(Path(plug_in_dir))
    click.echo("run_protocol")
    click.echo(plugin)

@click.command()
@click.option('--output', default="protocol_spec.txt", help="Name of specification output file", show_default=True)
@click.argument('plug_in_dir')
@click.argument('data_dir')
@click.argument('runsheet_path')
def spec(plug_in_dir, output, data_dir, runsheet_path):
    plugin = load_plugin(Path(plug_in_dir))
    output = Path(output)
    data_dir = Path(data_dir)
    runsheet_path = Path(runsheet_path)
    click.echo(f"Generating specification of validation protocol and outputting to file: '{output}'")\
    
    datasystem = load_data(
        config=plugin.config,
        root_path=data_dir,
        runsheet_path=runsheet_path,
    )

    vp = plugin.protocol.validate(
        datasystem.dataset,
        report_args={"include_skipped": True},
        defer_run=True,
    )

    specification = vp.queued_checks(
                long_description = True,
                CHECK_PREFIX = "",
                COMPONENT_PREFIX = " ",
                INDENT_CHAR = "#",
                WRAP_COMPONENT_NAME_CHAR = "",
                include_checks_counters = False
                )

    click.echo(specification)

    with open(output, "w") as f:
        f.write(specification)

    click.echo(f"Saved specification to {output}")

validation.add_command(run)
validation.add_command(spec)


from pathlib import Path
import re
import itertools

import click
from loguru import logger
import yaml

from dp_tools.glds_api import commons, isa
from dp_tools.core.files import isa_archive

@click.group()
def data_assets():
    pass

def product_dict(**kwargs):
    keys = kwargs.keys()
    vals = kwargs.values()
    for instance in itertools.product(*vals):
        yield dict(zip(keys, instance))

def matches_template(query_path: Path, template_paths: list[str], template_config: dict, is_directory: bool):
    logger.trace(query_path)
    logger.trace(template_paths)
    logger.trace(template_config)
    logger.trace(is_directory)

    template_str = str(Path(template_paths[0]).joinpath(*[Path(p) for p in template_paths[1:]]))

    logger.trace(f"template_str: {template_str}")


    # See if any template exist
    if not re.search('{.*}', template_str):
        logger.trace("No template sections, using exact filename")
        possible_filenames = [template_str]
    else:
        # Create possible filenames by filling in template
        possible_filenames = [template_str.format(**template_value_combo) for template_value_combo in product_dict(**template_config)]
        logger.trace(possible_filenames )

    for possible_filename in possible_filenames:
        # Guard clause: Directory style assets
        if query_path.is_relative_to(Path(possible_filename)) and is_directory:
            return True

        if query_path.match(possible_filename):
            return True

    #logger.trace(f"Matched file: {path_file} to template key: '{key}' with local path template: '{template_path}'")
    return False


def load_config(yaml_file: str):
    logger.info(f"Found existing yaml file: loading {yaml_file}")
    with open(yaml_file, "r") as f:
        config = yaml.safe_load(f)
    return config

@click.command()
@click.argument("root-dir")
@click.argument("template-values-yaml")
@click.argument("exclude-patterns", nargs=-1)
@click.option("--data-asset-yaml", default="data_assets.yaml")
def generate_config(root_dir: str, exclude_patterns: list[str], template_values_yaml:str,  data_asset_yaml: str):
    if Path(data_asset_yaml).exists():
        config = load_config(data_asset_yaml)
    else:
        logger.info(f"Generating new yaml file at: {data_asset_yaml}")
        config = {
            "data assets": {}
        }
        
    data_assets = config["data assets"]

    template_values = load_config(template_values_yaml)

    logger.info(f"Loaded existing configuration file with '{len(data_assets)}' data assets")
    logger.trace(f"Full data assets loaded: {data_assets}")
    path_root_dir = Path(root_dir)

    logger.info("Generate config starting..")

    logger.info(f"Reading files from {path_root_dir}")

    # Get all file paths as relative to path_root_dir
    path_file_all = [f.relative_to(path_root_dir) for f in path_root_dir.rglob("*") if f.is_file()] # list to evaluate generator

    # Apply any requested filters
    path_file_filtered: list[Path] = path_file_all.copy()
    to_remove: set[Path] = set()
    pattern_hit_counter: dict[str, int] = dict()
    for pattern in exclude_patterns:
        pattern_hit_counter[pattern] = 0
        logger.trace(f"Checking pattern: {pattern}")
        for f in path_file_filtered:
            if f.match(pattern):
                pattern_hit_counter[pattern] = pattern_hit_counter[pattern] + 1
                logger.trace(f"Filtered out {f} due to matching excluded pattern: {pattern}")
                to_remove.add(f)
                continue
    
    # Note: we remove after checking all files.
    #   Removal from list during iteration causes an iteration bug
    for f in to_remove:
        path_file_filtered.remove(f)


    # Logging re: filtering results
    for pattern, filter_count in pattern_hit_counter.items():
        if filter_count != 0:
            logger.debug(f"Filtered out {filter_count} files that matched this pattern: '{pattern}'")
        else:
            logger.warning(f"Filtered out 0 files based on this pattern: '{pattern}'. Double check the exclusion pattern is correct.")


    if len(path_file_all) != len(path_file_filtered):
        logger.info(f"Filtered out a total of {len(path_file_all) - len(path_file_filtered)} of {len(path_file_all)} to a total of {len(path_file_filtered)} files across all exclude patterns: {exclude_patterns}")

    # Trace logging for all files that were loaded
    def template_path_to_glob_path(s: str) -> str:
        """Convert template style string to glob key

        :param s: A template string. E.g. {{sample}}_results.out
        :type s: str
        :return: A equivalent glob string. E.g. *_results.out
        :rtype: str
        """
        glob_path = re.sub("{{.*?}}","*",s)
        logger.trace(f"Converted {s} to {glob_path}")
        return glob_path
    path_file_unassigned: list[Path] = path_file_filtered.copy()
    logger.info(f"Checking {len(path_file_unassigned)} for matching data asset specifications")

    while path_file_unassigned:
        for i, path_file in enumerate(path_file_unassigned):
            path_has_match = False
            logger.debug(f"Checking if {path_file} matches existing keys. File Number {i+1} of {len(path_file_unassigned)}")
            for key, meta_asset in data_assets.items():
                logger.trace(meta_asset)
                if matches_template(path_file, meta_asset['local location'], template_values, meta_asset.get("is directory")):
                    path_file_unassigned.remove(path_file)
                    path_has_match = True
                    continue

            if path_has_match:
                continue
            
            new_key = ""
            while True:
                if new_key == "":
                    new_key = click.prompt(f"No data asset key found for {path_file} {i+1} of {len(path_file_unassigned)}\n\tPlease enter data asset key")
                    if new_key in data_assets:
                        click.echo(f"This data key name is already used! Must supply a unique asset key name.")
                        continue
                new_template_path = click.prompt(f"For key: {new_key} target file: {path_file}\n\tPlease enter new local path template")
                treat_as_directory = click.confirm("Is this a directory data asset?")
                if not matches_template(path_file, new_template_path.split('/'), template_values, treat_as_directory):
                    click.echo(f"'{new_template_path}' does not match {path_file}!")
                    continue
                else:
                    data_assets[new_key] = {
                            "local location": new_template_path.split("/"),
                            "is directory": treat_as_directory
                        }
                    # Save to file each new data asset added
                    with open(data_asset_yaml, "w") as f:
                        yaml.dump(config, f)
                    logger.success(f"Added {new_key} to configuration!")
                    break # Get out of new_key loop
                    

    
    


data_assets.add_command(generate_config)


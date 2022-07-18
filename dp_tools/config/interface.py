""" Functions that parse configuration files """
from typing import Union
from pathlib import Path
import os
import yaml
import logging
log = logging.getLogger(__name__)
from importlib.metadata import files

ConfigVersion = tuple[str, str]
""" Denotes a specific prepackaged configuration file with (<assay>,<version>), e.g. ("bulkRNAseq","Latest") """

ConfigSelection = Union[ConfigVersion, Path]
""" Specifies a configuration file, either prepackaged or by direct local path """

def load_config(config: ConfigSelection) -> dict:
    """Load yaml configuration file. Allows loading from either:
      - A prepackaged configuration file using a tuple of ('config_type','config_version') (e.g. ('bulkRNASeq','Latest'), ('microarray','0'))
      - A configuration file supplied as a Path object

    :param config: Configuration file to load
    :type config: Union[tuple[str,str], Path]
    :return: A dictionary of the full configuration
    :rtype: dict
    """
    match config:
        case tuple():
            conf_type, conf_version = config
            query_config_fn = f"{conf_type}_v{conf_version}.yaml"
            [resolved_config_path] = (p for p in files('dp_tools') if p.name == query_config_fn)
            log.info(f"Loading config (relative to package): {resolved_config_path}")
            with open(resolved_config_path.locate(), "r") as f:
                conf_full = yaml.safe_load(f)
        case Path():
            log.info(f"Loading config (direct path): {config}")
            conf_full = yaml.safe_load(config.open())

    log.debug(f"Final config loaded: {conf_full}")

    return conf_full


def get_data_asset_template(key: str, config: ConfigSelection) -> str:
    """Fetches the canonical processing output location for a given key

    :param key: Data asset section key name
    :type key: str
    :return: A string denoting a path to the data asset. It may contain unfilled template sections which can be filled with the 'format' method
    :rtype: str
    """
    ...

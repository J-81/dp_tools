import os
from pathlib import Path
from typing import Union
import pkg_resources
from warnings import warn

import yaml
import logging

log = logging.getLogger(__name__)


def load_full_config(config: Union[str, Path]) -> dict:
    warn(
        "Calls to this function should be migrated to 'load_config'; version=2.0.0",
        DeprecationWarning,
        stacklevel=2,
    )
    if isinstance(config, str):
        resolved_config_path = os.path.join(
            "..", "config", f"bulkRNASeq_v{config}.yaml"
        )
        log.info(f"Loading full config (relative to package): {resolved_config_path}")
        conf_full = yaml.safe_load(
            pkg_resources.resource_string(__name__, resolved_config_path)
        )
    elif isinstance(config, Path):
        log.info(f"Loading config (direct path): {config}")
        conf_full = yaml.safe_load(config.open())

    # validate with schema
    # config_schema = Schema()

    log.debug(f"Final config loaded: {conf_full}")

    return conf_full

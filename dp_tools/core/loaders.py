from pathlib import Path
import logging

log = logging.getLogger(__name__)

from dp_tools.core.configuration import load_config
from dp_tools.core.entity_model2 import DataSystem, dataSystem_from_runsheet


def load_data(
    config: tuple[str, str] | Path,
    root_path: Path,
    runsheet_path: Path,
    key_sets: tuple[str] = None,
    keys: list[str] = None,
) -> DataSystem:
    # Load config
    conf = load_config(config)

    # determine keys
    if keys is None:
        keys = list()

    if key_sets is not None:
        for key_set in key_sets:
            keys.extend(conf["data asset sets"][key_set])
    log.info(f"Attempting to load data for {len(keys)} data asset keys")
    log.debug(f"Attempting to load data for these data asset keys: {keys}")

    # create data system
    dataSystem = dataSystem_from_runsheet(runsheet_path=runsheet_path)

    # create dataset
    dataset = dataSystem.dataset_from_runsheet(runsheet_path)

    # Load data assets configuration
    conf_data_assets = conf["data assets"]

    # Check all requested data assets keys are present
    assert set(keys).issubset(
        set(conf_data_assets)
    ), f"Could not find {set(keys) - set(conf_data_assets)} in data asset keys"

    # Load Metadata from runsheet
    for key in keys:
        dataset.load_data_asset(
            data_asset_config=conf_data_assets[key], name=key, root_dir=root_path
        )

    return dataSystem

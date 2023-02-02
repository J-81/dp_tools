import logging

from dp_tools.core.loaders import load_data
from dp_tools.core.configuration import available_data_asset_keys

import pytest


def test_list_keys_in_config():
    keys = available_data_asset_keys(config=("bulkRNASeq", "1"))
    assert len(keys) == 74


def test_key_based_loader_GLDS48(root_test_dir):

    ds = load_data(
        key_sets=["is single end full"],
        config=("bulkRNASeq", "Latest"),
        root_path=(root_test_dir / "GLDS-48"),
        runsheet_path=(
            root_test_dir / "GLDS-48/Metadata/GLDS-48_bulkRNASeq_v1_runsheet.csv"
        ),
    )


def test_key_based_loader_GLDS194(root_test_dir, caplog):
    caplog.set_level(logging.INFO)

    load_data(
        key_sets=["is paired end full", "has ercc"],
        config=("bulkRNASeq", "Latest"),
        root_path=(root_test_dir / "GLDS-194"),
        runsheet_path=(
            root_test_dir / "GLDS-194/Metadata/GLDS-194_bulkRNASeq_v1_runsheet.csv"
        ),
    )

def test_key_based_loader_bad_keys(root_test_dir):
    with pytest.raises(AssertionError):
        load_data(
            keys=["typoassetkey"],
            config=("bulkRNASeq", "Latest"),
            root_path=(root_test_dir / "GLDS-194"),
            runsheet_path=(
                root_test_dir / "GLDS-194/Metadata/GLDS-194_bulkRNASeq_v1_runsheet.csv"
            ),
        )

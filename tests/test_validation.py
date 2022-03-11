""" Tests for validation report results, relies on test for loaders passing """
from pathlib import Path
import os

from dp_tools.bulkRNASeq.loaders import load_from_bulk_rnaseq_raw_dir
from dp_tools.core.entity_model import Flag

# set for testing
TEST_DIR = Path(os.environ["TEST_ASSETS_DIR"])


def test_bulk_rnaseq_raw_data_validation(caplog):
    """ This tests validation as it would be run on dataset after demultiplexing """

    target_data_dir = TEST_DIR / "GLDS-194"

    caplog.set_level(0)
    ds = load_from_bulk_rnaseq_raw_dir(
        target_data_dir
    )

    ds.validate()

    assert isinstance(ds.validation_report['dataset'], dict)
    assert isinstance(ds.validation_report['sample'], dict)
    assert isinstance(ds.validation_report['component'], dict)

    assert [isinstance(flag, Flag) for flag in ds.validation_report['dataset'].values()]
    assert [isinstance(flag, Flag) for flag in ds.validation_report['sample'].values()]
    assert [isinstance(flag, Flag) for flag in ds.validation_report['component'].values()]
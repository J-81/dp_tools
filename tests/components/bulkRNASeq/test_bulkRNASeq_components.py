

import os
from pathlib import Path
from dp_tools.bulkRNASeq.loaders import load_BulkRNASeq_STAGE_00
import pytest

# set for testing
TEST_DIR = Path(os.environ["TEST_ASSETS_DIR"])

@pytest.fixture
def dataSystem():
    target_data_dir = TEST_DIR / "GLDS-207_PostDemultiplex"

    ds = load_BulkRNASeq_STAGE_00(
        target_data_dir,
        dataSystem_name = "GLDS-207"
    )

    return ds


def test_rawReadsComponent(dataSystem):
    # iterate through all raw read components
    for sample in dataSystem.dataset.samples.values():
        rrComp = sample.rawForwardReads

        assert rrComp.count == 100
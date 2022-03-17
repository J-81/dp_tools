import os
from pathlib import Path
from dp_tools.bulkRNASeq.loaders import load_BulkRNASeq_STAGE_00
import pytest


def test_rawReadsComponent(glds194_dataSystem_STAGE00):
    # iterate through all raw read components
    for sample in glds194_dataSystem_STAGE00.dataset.samples.values():

        assert sample.rawForwardReads.total_sequences == 200
        assert sample.rawReverseReads.total_sequences == 200
        assert sample.rawForwardReads.total_sequences == 200
        assert sample.rawReverseReads.total_sequences == 200


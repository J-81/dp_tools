""" Functions for loading into data model objects 

Loaders SHOULD:
  - run validation on the dataSystem object
  - return a dataSystem object
"""
from pathlib import Path
import logging
from typing import List, Protocol, Union

import pandas as pd

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)

from dp_tools.model import (
    BaseComponent,
    BaseSample,
    BaseDataSystem,
    BaseDataset,
    BulkRNASeqDataset,
    BulkRNASeqSample,
    DataFile,
    GLDSDataSystem,
    ReadsComponent,
)
from dp_tools.metadata import BulkRNASeqMetadata
from dp_tools.locaters import RawFastq


def load_from_bulk_rnaseq_raw_dir(root_path: Path, metadata: BulkRNASeqMetadata):
    # assert protocols are implemented
    assert isinstance(
        metadata, BulkRNASeqMetadata
    ), f"{metadata} does not implement proper protcol for {BulkRNASeqMetadata}"

    # create datasystem
    log.info(f"Attempting to load data model for raw directory: {str(root_path)}")
    dataSystem = GLDSDataSystem( base=BaseDataSystem(name=root_path.name))
    # initate locaters on the root path
    rawFastq = RawFastq(search_root=root_path)

    # create dataset
    dataset = BulkRNASeqDataset(base=BaseDataset(name=root_path.name))

    dataSystem.attach_dataset(dataset)

    # create samples
    for sample_name in metadata.samples():
        sample = BulkRNASeqSample(BaseSample(name=sample_name))
        if metadata.paired_end:
            raw_fwd_reads = ReadsComponent(
                base=BaseComponent(description="Raw Reverse Reads"),
                fastqGZ=DataFile(rawFastq.find(sample_name, type=("Raw", "Forward"))),
            )
            raw_rev_reads = ReadsComponent(
                base=BaseComponent(description="Raw Reverse Reads"),
                fastqGZ=DataFile(
                    path=rawFastq.find(sample_name, type=("Raw", "Reverse"))
                ),
            )
            sample.attach_component(raw_fwd_reads, attr="rawForwardReads")
            sample.attach_component(raw_rev_reads, attr="rawReverseReads")
        else:
            raw_reads = ReadsComponent(
                base=BaseComponent(description="Raw Reads"),
                fastqGZ=DataFile(
                    path=rawFastq.find(sample_name, type=("Raw", "Forward"))
                ),
            )
            sample.attach_component(raw_reads, attr="rawReads")
        # attach components
        dataset.attach_sample(sample)

    # run validators
    dataSystem.validate()

    return dataSystem

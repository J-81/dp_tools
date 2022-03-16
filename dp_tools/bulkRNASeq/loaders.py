""" Functions for loading into data model objects 

Loaders SHOULD:
  - run validation on the dataSystem object
  - return a dataSystem object
"""
import os
from pathlib import Path
import logging
from typing import List, Protocol, Union

import pandas as pd

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)

from dp_tools.core.entity_model import (
    BaseComponent,
    BaseSample,
    BaseDataSystem,
    BaseDataset,
    DataDir,
    DataFile,
    GLDSDataSystem,
)

from dp_tools.bulkRNASeq.entity import BulkRNASeqDataset, BulkRNASeqSample
from dp_tools.bulkRNASeq.locaters import FastqcReport, MultiQCDir, RawFastq, Runsheet
from dp_tools.components import RawReadsComponent, BulkRNASeqMetadataComponent


def load_from_bulk_rnaseq_raw_dir(root_path: Path, dataSystem_name: str = None):
    # ensure root path exists!
    if not root_path.is_dir():
        raise FileNotFoundError(f"Root path doesn't exist!: {root_path}")

    # if not specified, assume root_path name is the dataSystem name
    if not dataSystem_name:
        dataSystem_name = root_path.name
    # create datasystem
    log.info(f"Attempting to load data model for raw directory: {str(root_path)}")
    dataSystem = GLDSDataSystem(base=BaseDataSystem(name=dataSystem_name))

    # initate locaters on the root path
    rawFastq = RawFastq(search_root=root_path)
    runsheet = Runsheet(search_root=root_path)
    readsMQC = MultiQCDir(search_root=root_path)
    fastqcReport = FastqcReport(search_root=root_path)

    # create dataset
    dataset = BulkRNASeqDataset(base=BaseDataset(name=dataSystem_name))

    dataSystem.attach_dataset(dataset)
    # attach dataset components
    dataSystem.dataset.attach_component(
        BulkRNASeqMetadataComponent(
            base=BaseComponent(description="Metadata in a runsheet csv file"),
            runsheet=DataFile(runsheet.find(datasystem_name=dataSystem_name)),
        ),
        attr="metadata",
    )

    # alias metadata for convenience
    metadata = dataSystem.dataset.all_components["metadata"]

    # create shared sample datafiles
    datf_readsMQC = DataDir(
        readsMQC.find(
            rel_dir=Path.joinpath(Path("00-RawData"), Path("FastQC_Reports")),
            mqc_label="raw",
        )
    )

    # create samples
    for sample_name in metadata.samples:
        sample = BulkRNASeqSample(BaseSample(name=sample_name))
        if metadata.paired_end:
            raw_fwd_reads = RawReadsComponent(
                base=BaseComponent(description="Raw Forward Reads"),
                fastqGZ=DataFile(rawFastq.find(sample_name, type=("Raw", "Forward"))),
                multiQCDir=datf_readsMQC,
                fastqcReportHTML=DataFile(path=fastqcReport.find(sample_name=sample_name, ext="html", type=("Raw", "Forward"))),
                fastqcReportZIP=DataFile(path=fastqcReport.find(sample_name=sample_name, ext="zip", type=("Raw", "Forward")))
            )
            raw_rev_reads = RawReadsComponent(
                base=BaseComponent(description="Raw Reverse Reads"),
                fastqGZ=DataFile(
                    path=rawFastq.find(sample_name, type=("Raw", "Reverse"))
                ),
                multiQCDir=datf_readsMQC,
                fastqcReportHTML=DataFile(path=fastqcReport.find(sample_name=sample_name, ext="html", type=("Raw", "Reverse"))),
                fastqcReportZIP=DataFile(path=fastqcReport.find(sample_name=sample_name, ext="zip", type=("Raw", "Reverse")))
            )
            sample.attach_component(raw_fwd_reads, attr="rawForwardReads")
            sample.attach_component(raw_rev_reads, attr="rawReverseReads")
        else:
            raw_reads = RawReadsComponent(
                base=BaseComponent(description="Raw Reads"),
                fastqGZ=DataFile(
                    path=rawFastq.find(sample_name, type=("Raw", "Forward"))
                ),
                multiQCDir=datf_readsMQC,
                fastqcReportHTML=DataFile(path=fastqcReport.find(sample_name=sample_name, ext="html", type=("Raw", "Forward"))),
                fastqcReportZIP=DataFile(path=fastqcReport.find(sample_name=sample_name, ext="zip", type=("Raw", "Forward")))
            )
            sample.attach_component(raw_reads, attr="rawReads")
        # attach components
        dataset.attach_sample(sample)

    return dataSystem

from collections import defaultdict
import enum
from typing import Dict, List
import logging
log = logging.getLogger(__name__)

from dp_tools.bulkRNASeq.entity import BulkRNASeqDataset, BulkRNASeqSample
from dp_tools.core.check_model import Flag, VVProtocol
from dp_tools.bulkRNASeq.checks import (
    COMPONENT_RAWREADS_0001,
    COMPONENT_TRIMREADS_0001,
    DATASET_RAWREADS_0001,
    DATASET_TRIMREADS_0001,
    SAMPLE_RAWREADS_0001,
)
from dp_tools.core.entity_model import TemplateComponent
from dp_tools.components import BulkRNASeqMetadataComponent, RawReadsComponent, TrimReadsComponent

class STAGE(enum.Enum):
    Demultiplexed = 0
    Reads_PreProcessed = 1
    GenomeAligned = 2
    RSeQCAnalysis = 2.01
    GeneCounted = 3
    DGE = 4

    # allow comparing stages
    # used to check if a sufficient stage has been achieved
    def __ge__(self, other):
        return self.value >= other.value


class BulkRNASeq_VVProtocol_RawData(VVProtocol):
    expected_dataset_class = BulkRNASeqDataset
    stage: STAGE = STAGE.Demultiplexed

    # init all checks
    dataset_rawreads_0001 = DATASET_RAWREADS_0001()
    sample_rawreads_0001 = SAMPLE_RAWREADS_0001()
    component_rawreads_0001 = COMPONENT_RAWREADS_0001()
    component_trimreads_0001 = COMPONENT_TRIMREADS_0001()

    def validate_dataset(self) -> Dict[BulkRNASeqDataset, List[Flag]]:
        flags: Dict[BulkRNASeqDataset, List[Flag]] = defaultdict(list)
        if self.stage >= STAGE.Demultiplexed: flags[self.dataset].append(self.dataset_rawreads_0001.validate(self.dataset))
        if self.stage >= STAGE.Reads_PreProcessed: flags[self.dataset].append(self.dataset_rawreads_0001.validate(self.dataset))
        return flags

    def validate_samples(self) -> Dict[BulkRNASeqSample, List[Flag]]:
        flags: Dict[BulkRNASeqSample, List[Flag]] = defaultdict(list)
        for sample in self.dataset.samples.values():
            if self.stage >= STAGE.Demultiplexed: flags[sample].append(self.sample_rawreads_0001.validate(sample=sample))
        return flags

    def validate_components(self) -> Dict[TemplateComponent, List[Flag]]:
        flags: Dict[TemplateComponent, List[Flag]] = defaultdict(list)
        # iterate through all components by level
        for component in self.dataset.all_non_empty_components_recursive:
            log.info(f"Validating component: {component} with type {type(component)}")
            match component:
                case RawReadsComponent():
                    if self.stage >= STAGE.Demultiplexed: flags[component].append(self.component_rawreads_0001.validate(component))
                case TrimReadsComponent():
                    if self.stage >= STAGE.Reads_PreProcessed: flags[component].append(self.component_trimreads_0001.validate(component))
                case BulkRNASeqMetadataComponent():
                    flags[component] = list()
                case _:
                    raise TypeError(f"Encountered unhandled component type in VV: {component} with type {type(component)}")
        return flags

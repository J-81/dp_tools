from collections import defaultdict
import enum
from typing import Dict, List
import logging

from dp_tools.components.components import DatasetGeneCounts, DifferentialGeneExpression, GeneCounts, GenomeAlignments, NormalizedGeneCounts, RSeQCAnalysis
log = logging.getLogger(__name__)

from dp_tools.bulkRNASeq.entity import BulkRNASeqDataset, BulkRNASeqSample
from dp_tools.core.check_model import Flag, VVProtocol
from dp_tools.bulkRNASeq.checks import (
    COMPONENT_GENOMEALIGNMENTS_0001,
    COMPONENT_RAWREADS_0001,
    COMPONENT_TRIMREADS_0001,
    DATASET_GENECOUNTS_0001,
    DATASET_GENOMEALIGNMENTS_0001,
    DATASET_METADATA_0001,
    DATASET_RAWREADS_0001,
    DATASET_RSEQCANALYSIS_0001,
    DATASET_TRIMREADS_0001,
    SAMPLE_RAWREADS_0001,
    SAMPLE_TRIMREADS_0001,
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


class BulkRNASeq_VVProtocol(VVProtocol):
    expected_dataset_class = BulkRNASeqDataset

    # init all checks
    # dataset level
    dataset_metadata_0001 = DATASET_METADATA_0001()
    dataset_rawReads_0001 = DATASET_RAWREADS_0001()
    dataset_trimReads_0001 = DATASET_TRIMREADS_0001()
    dataset_genomeAlignments_0001 = DATASET_GENOMEALIGNMENTS_0001()
    dataset_rseqcAnalysis_0001 = DATASET_RSEQCANALYSIS_0001()
    dataset_geneCounts = DATASET_GENECOUNTS_0001()
    
    # sample level
    sample_rawReads_0001 = SAMPLE_RAWREADS_0001()
    sample_trimReads_0001 = SAMPLE_TRIMREADS_0001()
    #sample_geneCounts = SAMPLE_GENECOUNTS_0001()

    # component level
    component_rawReads_0001 = COMPONENT_RAWREADS_0001()
    component_trimReads_0001 = COMPONENT_TRIMREADS_0001()
    component_genomeAlignments_0001 = COMPONENT_GENOMEALIGNMENTS_0001()

    def __init__(self, dataset: BulkRNASeqDataset, stage: STAGE):
        super().__init__(dataset)
        self.stage = stage

    def validate_dataset(self) -> Dict[BulkRNASeqDataset, List[Flag]]:
        flags: Dict[BulkRNASeqDataset, List[Flag]] = defaultdict(list)
        if self.stage >= STAGE.Demultiplexed: flags[self.dataset].append(self.dataset_metadata_0001.validate(self.dataset))
        if self.stage >= STAGE.Demultiplexed: flags[self.dataset].append(self.dataset_rawReads_0001.validate(self.dataset))
        if self.stage >= STAGE.Reads_PreProcessed: flags[self.dataset].append(self.dataset_trimReads_0001.validate(self.dataset))
        if self.stage >= STAGE.GenomeAligned: flags[self.dataset].append(self.dataset_genomeAlignments_0001.validate(self.dataset))
        if self.stage >= STAGE.RSeQCAnalysis: flags[self.dataset].append(self.dataset_rseqcAnalysis_0001.validate(self.dataset))
        if self.stage >= STAGE.GeneCounted: flags[self.dataset].append(self.dataset_geneCounts.validate(self.dataset))
        return flags

    def validate_samples(self) -> Dict[BulkRNASeqSample, List[Flag]]:
        flags: Dict[BulkRNASeqSample, List[Flag]] = defaultdict(list)
        for sample in self.dataset.samples.values():
            if self.stage >= STAGE.Demultiplexed: flags[sample].append(self.sample_rawReads_0001.validate(sample=sample))
            if self.stage >= STAGE.Reads_PreProcessed: flags[sample].append(self.sample_trimReads_0001.validate(sample=sample))
        return flags

    def validate_components(self) -> Dict[TemplateComponent, List[Flag]]:
        flags: Dict[TemplateComponent, List[Flag]] = defaultdict(list)
        # iterate through all components by level
        for component in self.dataset.all_non_empty_components_recursive:
            log.info(f"Validating component: {component} with type {type(component)}")
            match component:
                case RawReadsComponent():
                    if self.stage >= STAGE.Demultiplexed: flags[component].append(self.component_rawReads_0001.validate(component))
                case TrimReadsComponent():
                    if self.stage >= STAGE.Reads_PreProcessed: flags[component].append(self.component_trimReads_0001.validate(component))
                case BulkRNASeqMetadataComponent():
                    flags[component] = list()
                case GenomeAlignments():
                    if self.stage >= STAGE.GenomeAligned: flags[component].append(self.component_genomeAlignments_0001.validate(component))
                case RSeQCAnalysis():
                    pass # no component level checks implemented
                case GeneCounts():
                    pass # no component level checks implemented
                case DatasetGeneCounts():
                    pass # no component level checks implemented
                case NormalizedGeneCounts():
                    pass # no component level checks implemented
                case DifferentialGeneExpression():
                    pass # no component level checks implemented
                case _:
                    raise TypeError(f"Encountered unhandled component type in VV: {component} with type {type(component)}")
        return flags

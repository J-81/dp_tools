from collections import defaultdict
from typing import Dict, List
import logging
log = logging.getLogger(__name__)

from dp_tools.bulkRNASeq.entity import BulkRNASeqDataset, BulkRNASeqSample
from dp_tools.core.check_model import Flag, VVProtocol
from dp_tools.bulkRNASeq.checks import (
    COMPONENT_RAWREADS_0001,
    COMPONENT_TRIMREADS_0001,
    SAMPLE_RAWREADS_0001,
)
from dp_tools.core.entity_model import TemplateComponent
from dp_tools.components import BulkRNASeqMetadataComponent, RawReadsComponent, TrimReadsComponent


class BulkRNASeq_VVProtocol_RawData(VVProtocol):
    expected_dataset_class = BulkRNASeqDataset

    def validate_dataset(self) -> Dict[BulkRNASeqDataset, List[Flag]]:
        flags: Dict[BulkRNASeqDataset, List[Flag]] = defaultdict(list)
        flags[self.dataset] = list()
        return flags

    def validate_samples(self) -> Dict[BulkRNASeqSample, List[Flag]]:
        flags: Dict[BulkRNASeqSample, List[Flag]] = defaultdict(list)
        for sample in self.dataset.samples.values():
            flags[sample].append(SAMPLE_RAWREADS_0001.validate(sample=sample))
        return flags

    def validate_components(self) -> Dict[TemplateComponent, List[Flag]]:
        flags: Dict[TemplateComponent, List[Flag]] = defaultdict(list)
        # iterate through all components
        for component in self.dataset.components:
            log.info(f"Validating component: {component} with type {type(component)}")
            match component:
                case RawReadsComponent():
                    flags[component].append(COMPONENT_RAWREADS_0001.validate(component))
                case TrimReadsComponent():
                    flags[component].append(COMPONENT_TRIMREADS_0001.validate(component))
                case BulkRNASeqMetadataComponent():
                    flags[component] = list()
                case _:
                    raise TypeError(f"Encountered unhandled component type in VV: {component} with type {type(component)}")
        return flags

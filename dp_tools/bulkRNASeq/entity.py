############################################################################################
# BULK RNASEQ SPECIFIC
############################################################################################
from dataclasses import dataclass, field
import enum
from typing import List, Union
from dp_tools.bulkRNASeq.checks import SAMPLE_RAWREADS_0001
from dp_tools.core.check_model import Flag

from dp_tools.core.entity_model import (
    BaseDataset,
    BaseSample,
    CanAttachComponents,
    EmptyComponent,
    TemplateDataset,
    TemplateSample,
)
from dp_tools.core.model_commons import strict_type_checks
from dp_tools.components import BulkRNASeqMetadataComponent, ReadsComponent


class ASSAY(enum.Enum):
    BulkRNASeq = 1


@dataclass(eq=False)
class BulkRNASeqSample(TemplateSample, CanAttachComponents):
    """ Abstract class for samples """

    # composition for all samples
    base: BaseSample = field(repr=False)
    assay_type: str = ASSAY.BulkRNASeq.name

    # used for paired end
    rawForwardReads: Union[EmptyComponent, ReadsComponent] = field(
        default_factory=EmptyComponent
    )
    rawReverseReads: Union[EmptyComponent, ReadsComponent] = field(
        default_factory=EmptyComponent
    )
    trimForwardReads: Union[EmptyComponent, ReadsComponent] = field(
        default_factory=EmptyComponent
    )
    trimReverseReads: Union[EmptyComponent, ReadsComponent] = field(
        default_factory=EmptyComponent
    )

    # used for single end
    rawReads: Union[EmptyComponent, ReadsComponent] = field(
        default_factory=EmptyComponent
    )
    TrimReads: Union[EmptyComponent, ReadsComponent] = field(
        default_factory=EmptyComponent
    )

    def __post_init__(self):
        pass

    def validate(self):
        flags = list()
        strict_type_checks(self)
        # additional checks advised
        flags.append(SAMPLE_RAWREADS_0001.validate(self))

        return flags  # may return empty list


@dataclass(eq=False)
class BulkRNASeqDataset(TemplateDataset, CanAttachComponents):

    base: BaseDataset = field(repr=False)
    assay_type: str = ASSAY.BulkRNASeq.name
    expected_sample_class = BulkRNASeqSample

    metadata: Union[EmptyComponent, BulkRNASeqMetadataComponent] = field(
        default_factory=EmptyComponent
    )

    def __post_init__(self):
        self.base.name = f"{self.base.name}__{self.assay_type}"

    def validate(self) -> List[Flag]:
        strict_type_checks(self, exceptions=["samples"])
        return list()  # may return empty list

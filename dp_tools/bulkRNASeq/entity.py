############################################################################################
# BULK RNASEQ SPECIFIC
############################################################################################
from dataclasses import dataclass, field
import enum
from typing import Union

from dp_tools.core.entity_model import BaseDataset, BaseSample, CanAttachComponents, EmptyComponent, TemplateDataset, TemplateSample
from dp_tools.core.model_commons import strict_type_checks
from dp_tools.components import BulkRNASeqMetadataComponent, ReadsComponent

class ASSAY(enum.Enum):
    BulkRNASeq = 1


@dataclass(eq = False)
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
        strict_type_checks(self)
        # additional checks advised


@dataclass(eq = False)
class BulkRNASeqDataset(TemplateDataset, CanAttachComponents):

    base: BaseDataset = field(repr=False)
    assay_type: str = ASSAY.BulkRNASeq.name
    expected_sample_class = BulkRNASeqSample

    metadata: Union[EmptyComponent, BulkRNASeqMetadataComponent] = field(
        default_factory=EmptyComponent
    )

    def __post_init__(self):
        self.base.name = f"{self.base.name}__{self.assay_type}"

    def validate(self):
        strict_type_checks(self, exceptions=["samples"])

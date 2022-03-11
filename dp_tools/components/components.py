from dataclasses import dataclass, field
from typing import List, Tuple, Union

import logging
log = logging.getLogger(__name__)

import pandas as pd

from dp_tools.core.model_commons import strict_type_checks
from dp_tools.core.entity_model import BaseComponent, CanAttachEntity, DataDir, DataFile, TemplateComponent
from dp_tools.core.check_model import Flag


@dataclass(eq=False)
class ReadsComponent(TemplateComponent, CanAttachEntity):

    fastqGZ: DataFile
    base: BaseComponent = field(repr=False)
    # id: str = field(default_factory=get_id)
    multiQCDir: Union[DataDir, None] = field(default=None)
    fastqcReportHTML: Union[DataFile, None] = field(default=None)
    fastqcReportZIP: Union[DataFile, None] = field(default=None)
    trimmingReportTXT: Union[DataFile, None] = field(default=None)

    def validate(self) -> List[Flag]:
        log.critical("Pending validation implemention")
        return list()

@dataclass(eq=False)
class BulkRNASeqMetadataComponent(TemplateComponent, CanAttachEntity):

    base: BaseComponent = field(repr=False)
    runsheet: Union[DataFile, None] = field(default=None)
    samples: List = field(init=False) # NOTE: List[str] is a more precise type hint; however, this breaks strict type checking: https://bugs.python.org/issue44529
    paired_end: bool = field(init=False)

    def __post_init__(self):
        if self.runsheet:
            # extract runsheet as dataframe
            self.df = pd.read_csv(self.runsheet.path)
            self.samples = list(self.df['sample_name']) # explicit conversion from pandas series to standard list
            self.paired_end = bool(self.df["paired_end"].unique()[0]) # explicit conversion from numpy bool to standard bool
        
        strict_type_checks(self)

    def validate(self) -> List[Tuple]:
        log.critical("Pending validation implemention")
        return list()
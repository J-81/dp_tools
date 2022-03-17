from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Union

import logging

log = logging.getLogger(__name__)

import pandas as pd
import multiqc

from dp_tools.core.model_commons import strict_type_checks
from dp_tools.core.entity_model import (
    BaseComponent,
    DataDir,
    DataFile,
    TemplateComponent,
)
from dp_tools.core.utilites.multiqc_tools import format_as_dataframe, get_general_stats

class ReadsComponent(TemplateComponent):
    _count: int
    # TODO: refine this type hint
    _mqc_data: dict

    @property
    def count(self):
        if not getattr(self,'_count'):
            self._count = self.mqc_data
        return self._count
    
    @property
    def mqc_general_stats(self) -> Dict:
        if not getattr(self,'_mqc_data'):
            self._mqc_data = get_general_stats(multiqc.run(analysis_dir = [self.fastqcReportZIP.path], no_report=True, no_data_dir=True, plots_interactive=True))
        return self._mqc_data

@dataclass(eq=False)
class RawReadsComponent(ReadsComponent):

    fastqGZ: DataFile
    base: BaseComponent = field(repr=False)
    # id: str = field(default_factory=get_id)
    multiQCDir: Union[DataDir, None] = field(default=None)
    fastqcReportHTML: Union[DataFile, None] = field(default=None)
    fastqcReportZIP: Union[DataFile, None] = field(default=None)


@dataclass(eq=False)
class TrimReadsComponent(ReadsComponent):

    fastqGZ: DataFile
    base: BaseComponent = field(repr=False)
    # id: str = field(default_factory=get_id)
    multiQCDir: Union[DataDir, None] = field(default=None)
    fastqcReportHTML: Union[DataFile, None] = field(default=None)
    fastqcReportZIP: Union[DataFile, None] = field(default=None)
    trimmingReportTXT: Union[DataFile, None] = field(default=None)


@dataclass(eq=False)
class BulkRNASeqMetadataComponent(TemplateComponent):

    base: BaseComponent = field(repr=False)
    runsheet: Union[DataFile, None] = field(default=None)
    samples: List = field(
        init=False
    )  # NOTE: List[str] is a more precise type hint; however, this breaks strict type checking: https://bugs.python.org/issue44529
    paired_end: bool = field(init=False)

    def __post_init__(self):
        if self.runsheet:
            # extract runsheet as dataframe
            self.df = pd.read_csv(self.runsheet.path)
            self.samples = list(
                self.df["sample_name"]
            )  # explicit conversion from pandas series to standard list
            self.paired_end = bool(
                self.df["paired_end"].unique()[0]
            )  # explicit conversion from numpy bool to standard bool

        strict_type_checks(self)

from dataclasses import dataclass, field
from typing import Dict, List, OrderedDict, Union

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

# type hint definitions
SampleLevelDict = OrderedDict[  # samples level
    "str",  # denotes sample name
    OrderedDict["str", float],  # key:value level, e.g. 'percent_gc': 52.0
]

GeneralStatsDict = Dict["str", SampleLevelDict]  # modules level  # denotes module name

GeneralStatsDictForComponent = OrderedDict[  # samples level
    "str",  # denotes sample name
    OrderedDict["str", float],  # key:value level, e.g. 'percent_gc': 52.0
]


class ReadsComponent(TemplateComponent):
    _total_sequences: int
    # TODO: refine this type hint
    _mqc_data: dict
    _mqc_modules_limit: List[str] = [
        "FastQC"
    ]  # optional but speeds up multiqc processing time slightly and can constrain modules if desired

    def extract_mqc(self):
        """Extracts multiQC data which exposes attributes for each general stat as well as plot data as separate attributes."""
        self.mqc_general_stats
        # self.mqc_plot_data # TODO: Implement

    @property
    def mqc_general_stats(self) -> Dict:
        if not getattr(self, "_mqc_general_stats", None):
            general_stats_dict: GeneralStatsDict = get_general_stats(
                multiqc.run(
                    analysis_dir=[self.fastqcReportZIP.path],
                    no_report=True,
                    no_data_dir=True,
                    plots_interactive=True,  # ensure data is robustly populated (otherwise flat plots result in missing extractable data)
                    module=[
                        module.lower() for module in self._mqc_modules_limit
                    ],  # module names here are always lowercase
                )
            )
            # assert only one sample
            # then remove sample layer (but keep module layer for recognizing source multiqc parser)
            general_stats_for_this_component: GeneralStatsDictForComponent = dict()  # type: ignore
            for module in general_stats_dict:
                # assert only one sample
                assert (
                    len(general_stats_dict[module]) == 1
                ), "One and only one sample expected"
                # extract the one key,value dict and remove the sample layer by assignment overwrite
                general_stats_for_this_component[module] = list(
                    general_stats_dict[module].values()
                )[0]

                # set individual attributes
                for key, value in general_stats_for_this_component[module].items():
                    # check if this clashes with anything else
                    if getattr(self, key, None):
                        new_key = f"{module}_{key}"
                        log.warning(
                            f"Attribute: '{key}' already exists.  Prepending module to disambiguate: '{new_key}'"
                        )
                        key = new_key
                    log.debug(
                        f"Extracting {key} from multiQC module {module} into {self}"
                    )
                    setattr(self, key, value)

            # set whole dict here
            self._mqc_general_stats = general_stats_for_this_component

        return self._mqc_general_stats


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

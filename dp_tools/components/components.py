from dataclasses import dataclass, field
from pathlib import Path
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


@dataclass(eq=False)
class RawReadsComponent(TemplateComponent):

    fastqGZ: DataFile
    base: BaseComponent = field(repr=False)
    # id: str = field(default_factory=get_id)
    multiQCDir: Union[DataDir, None] = field(default=None)
    fastqcReportHTML: Union[DataFile, None] = field(default=None)
    fastqcReportZIP: Union[DataFile, None] = field(
        default=None, metadata={"mqc_parse": ["FastQC"]}
    )


@dataclass(eq=False)
class TrimReadsComponent(TemplateComponent):

    fastqGZ: DataFile
    base: BaseComponent = field(repr=False)
    # id: str = field(default_factory=get_id)
    multiQCDir: Union[DataDir, None] = field(default=None)
    fastqcReportHTML: Union[DataFile, None] = field(default=None)
    fastqcReportZIP: Union[DataFile, None] = field(
        default=None, metadata={"mqc_parse": ["FastQC"]}
    )
    # note: uses Cutadapt as Trim-Galore wraps cutadapt and Trim-Galore is not currently a multiQC module
    trimmingReportTXT: Union[DataFile, None] = field(
        default=None, metadata={"mqc_parse": ["Cutadapt"]}
    )


@dataclass(eq=False)
class BulkRNASeqMetadataComponent(TemplateComponent):

    base: BaseComponent = field(repr=False)
    runsheet: Union[DataFile, None] = field(default=None)
    samples: List = field(
        init=False
    )  # NOTE: List[str] is a more precise type hint; however, this breaks strict type checking: https://bugs.python.org/issue44529
    paired_end: bool = field(init=False)
    has_ercc: bool = field(init=False)

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
            self.has_ercc = bool(
                self.df["has_ERCC"].unique()[0]
            )  # explicit conversion from numpy bool to standard bool

        strict_type_checks(self)


@dataclass(eq=False)
class GenomeAlignments(TemplateComponent):

    base: BaseComponent = field(repr=False)
    alignedToTranscriptomeBam: Union[DataFile, None] = field(default=None)
    alignedSortedByCoordBam: Union[DataFile, None] = field(default=None)
    alignedSortedByCoordResortedBam: Union[DataFile, None] = field(default=None)
    alignedSortedByCoordResortedBamIndex: Union[DataFile, None] = field(default=None)
    logFinal: Union[DataFile, None] = field(
        default=None, metadata={"mqc_parse": ["STAR"]}
    )
    logProgress: Union[DataFile, None] = field(default=None)
    logFull: Union[DataFile, None] = field(default=None)
    sjTab: Union[DataFile, None] = field(default=None)
    multiQCDir: Union[DataDir, None] = field(default=None)


@dataclass(eq=False)
class GeneCounts(TemplateComponent):
    """ Based on RSEM output gene counts"""

    base: BaseComponent = field(repr=False)
    multiQCDir: Union[DataDir, None] = field(default=None)
    genesResults: Union[DataFile, None] = field(default=None)
    isoformsResults: Union[DataFile, None] = field(default=None)
    statDir: Union[DataDir, None] = field(
        default=None, metadata={"mqc_parse": ["Rsem"]}
    )


@dataclass(eq=False)
class DatasetGeneCounts(TemplateComponent):

    base: BaseComponent = field(repr=False)
    numNonZero: Union[DataFile, None] = field(default=None)
    unnormalizedCounts: Union[DataFile, None] = field(default=None)


@dataclass(eq=False)
class RSeQCAnalysis(TemplateComponent):

    base: BaseComponent = field(repr=False)
    geneBodyCoverageMultiQCDir: Union[DataDir, None] = field(default=None)
    geneBodyCoverageOut: Union[DataDir, None] = field(
        default=None, metadata={"mqc_parse": ["RSeQC"]}
    )
    inferExperimentMultiQCDir: Union[DataDir, None] = field(default=None)
    inferExperimentOut: Union[DataFile, None] = field(
        default=None, metadata={"mqc_parse": ["RSeQC"]}
    )
    innerDistanceMultiQCDir: Union[DataDir, None] = field(default=None)
    innerDistanceOut: Union[DataDir, None] = field(
        default=None, metadata={"mqc_parse": ["RSeQC"]}
    )
    readDistributionMultiQCDir: Union[DataDir, None] = field(default=None)
    readDistributionOut: Union[DataFile, None] = field(
        default=None, metadata={"mqc_parse": ["RSeQC"]}
    )


@dataclass(eq=False)
class NormalizedGeneCounts(TemplateComponent):

    base: BaseComponent = field(repr=False)
    erccNormalizedCountsCSV: Union[DataFile, None] = field(default=None)
    normalizedCountsCSV: Union[DataFile, None] = field(default=None)
    sampleTableCSV: Union[DataFile, None] = field(default=None)
    unnormalizedCountsCSV: Union[DataFile, None] = field(default=None)


@dataclass(eq=False)
class DifferentialGeneExpression(TemplateComponent):

    base: BaseComponent = field(repr=False)
    contrastsCSV: Union[DataFile, None] = field(default=None)
    annotatedTableCSV: Union[DataFile, None] = field(default=None)
    visualizationTableCSV: Union[DataFile, None] = field(default=None)
    visualizationPCATablCSV: Union[DataFile, None] = field(default=None)

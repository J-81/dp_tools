from dataclasses import dataclass, field
import itertools
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
    fastQCmultiQCDirZIP: Union[DataFile, None] = field(default=None)
    fastqcReportHTML: Union[DataFile, None] = field(default=None)
    fastqcReportZIP: Union[DataFile, None] = field(
        default=None, metadata={"mqc_parse": ["FastQC"]}
    )


@dataclass(eq=False)
class TrimReadsComponent(TemplateComponent):

    fastqGZ: DataFile
    base: BaseComponent = field(repr=False)
    # id: str = field(default_factory=get_id)
    fastQCmultiQCDirZIP: Union[DataFile, None] = field(default=None)
    fastqcReportHTML: Union[DataFile, None] = field(default=None)
    fastqcReportZIP: Union[DataFile, None] = field(
        default=None, metadata={"mqc_parse": ["FastQC"]}
    )
    # note: uses Cutadapt as Trim-Galore wraps cutadapt and Trim-Galore is not currently a multiQC module
    trimmingReportTXT: Union[DataFile, None] = field(
        default=None, metadata={"mqc_parse": ["Cutadapt"]}
    )
    trimmingMultiQCDirZIP: Union[DataFile, None] = field(default=None)


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

    @property
    def factor_groups(self) -> Dict[str, str]:
        """Returns sample wise dictionary of sample name to concatenated factor group.

        The concatenated group will be consistent with the BulkRNASeq Deseq2 script formatting for constrasts.
        This means:
         - replacing within factor 'spaces' with periods
         - wrapping the factor in quotes
         - joining factors into factor groups with ' & '
         - finally joining factors groups to derive pairwise-contrasts with 'v'

        """
        if getattr(self, "_factor_groups", None) == None:  # uncached
            if self.runsheet:
                # read as dataframe
                df = pd.read_csv(self.runsheet.path).set_index("sample_name")

                # find factors columns
                factor_cols = [c for c in df.columns if c.startswith("Factor Value[")]

                # create concatenated group column
                df["FactorGroup"] = (
                    "(" + df[factor_cols].agg(" & ".join, axis="columns") + ")"
                )

                self._factor_groups = df["FactorGroup"].to_dict()

            else:
                raise ValueError(
                    "No source of experimental factors available. Supported sources: 'runsheet'"
                )
        return self._factor_groups

    @property
    def contrasts(self) -> pd.DataFrame:
        """Returns sample wise dictionary of sample name to concatenated factor group.

        Creates a constrast dataframe in the same fashion as the DESEQ2 contrasts table:
        - determine all permutations of existing factor groups
        - format with header for each permutation in format: {factorGroup1}v{factorGroup2}
        - format with two rows for each column in format:
          - r1: factorGroup1.replace(' & ','...').replace(' ','.').replace('(',"").replace(')',"")
          - r2: factorGroup1.replace(' & ','...').replace(' ','.').replace('(',"").replace(')',"")
        """
        if getattr(self, "_contrasts_cached", None) == None:  # uncached
            # determine all permutations of existing factor groups
            permutations = list(
                itertools.permutations(set(self.factor_groups.values()), 2)
            )

            # format with header for each permutation
            headers = ["v".join(permute) for permute in permutations]

            # format with two rows for each column
            row1s = [
                permute[0]
                .replace(" & ", "...")
                .replace(" ", ".")
                .replace("(", "")
                .replace(")", "")
                for permute in permutations
            ]
            row2s = [
                permute[1]
                .replace(" & ", "...")
                .replace(" ", ".")
                .replace("(", "")
                .replace(")", "")
                for permute in permutations
            ]

            # format as dataframe
            self._contrasts = pd.DataFrame(
                {
                    header: [row1, row2]
                    for header, row1, row2 in zip(headers, row1s, row2s)
                }
            )
            # use this for cache check as:
            # 'ValueError: The truth value of a DataFrame is ambiguous. Use a.empty, a.bool(), a.item(), a.any() or a.all().'
            self._contrasts_cached = True

        return self._contrasts


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
    multiQCDirZIP: Union[DataFile, None] = field(default=None)


@dataclass(eq=False)
class GeneCounts(TemplateComponent):
    """Based on RSEM output gene counts"""

    base: BaseComponent = field(repr=False)
    multiQCDirZIP: Union[DataFile, None] = field(default=None)
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
    geneBodyCoverageMultiQCDirZIP: Union[DataFile, None] = field(default=None)
    geneBodyCoverageOut: Union[DataDir, None] = field(
        default=None, metadata={"mqc_parse": ["RSeQC"]}
    )
    inferExperimentMultiQCDirZIP: Union[DataFile, None] = field(default=None)
    inferExperimentOut: Union[DataFile, None] = field(
        default=None, metadata={"mqc_parse": ["RSeQC"]}
    )
    innerDistanceMultiQCDirZIP: Union[DataFile, None] = field(default=None)
    innerDistanceOut: Union[DataDir, None] = field(
        default=None, metadata={"mqc_parse": ["RSeQC"]}
    )
    readDistributionMultiQCDirZIP: Union[DataFile, None] = field(default=None)
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
    visualizationPCATableCSV: Union[DataFile, None] = field(default=None)

from dataclasses import dataclass, field
import itertools
from pathlib import Path
import tempfile
from typing import ClassVar, Dict, List, OrderedDict, Set, Union

import logging
import zipfile

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
    ISAarchive: Union[DataFile, None] = field(default=None)

    _ISA_INVESTIGATION_HEADERS: ClassVar[Set[str]] = {
        "ONTOLOGY SOURCE REFERENCE",
        "INVESTIGATION",
        "INVESTIGATION PUBLICATIONS",
        "INVESTIGATION CONTACTS",
        "STUDY",
        "STUDY DESIGN DESCRIPTORS",
        "STUDY PUBLICATIONS",
        "STUDY FACTORS",
        "STUDY ASSAYS",
        "STUDY PROTOCOLS",
        "STUDY CONTACTS",
    }

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
            [self.paired_end] = self.df["paired_end"].unique()
            self.paired_end = bool(self.paired_end)  # explicit conversion from numpy bool to standard bool
            [self.has_ercc] = self.df["has_ERCC"].unique()
            self.has_ercc = bool(self.has_ercc)  # explicit conversion from numpy bool to standard bool

        strict_type_checks(self, exceptions=["_ISA_INVESTIGATION_HEADERS"])

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

    def fetch_isa_files(self) -> Set[Path]:
        """Unzips the ISA archive in a temporary directory and reports files

        :return: List of individual file paths
        :rtype: List[Path]
        """
        if getattr(self, "_isa_files", None):
            self._isa_files: Set[Path]
            return self._isa_files

        assert self.ISAarchive, "No ISA archive data asset attached."

        temp_dir = tempfile.mkdtemp()
        log.debug(f"Extracting ISA Archive to temp directory: {temp_dir}")
        with zipfile.ZipFile(self.ISAarchive.path, "r") as zip_ref:
            zip_ref.extractall(temp_dir)

        self._isa_files = {f for f in Path(temp_dir).rglob("*") if f.is_file()}
        return self._isa_files

    def get_assays(self) -> Dict[str, Path]:
        """From the ISA Investigation file, extract assays mapped to tables

        :return: A mapping of assay name to table paths
        :rtype: Dict[str, Path]
        """

    @property
    def isa_investigation_subtables(self) -> dict[str, pd.DataFrame]:
        if cached := getattr(self, "_isa_investigation_subtables", None):
            return cached

        tables: dict[str, pd.DataFrame] = dict()

        # track sub table lines
        table_lines: List[list] = list()
        key: str = None  # type: ignore

        [i_file] = (f for f in self.fetch_isa_files() if f.name.startswith("i_"))
        with open(i_file, "r") as f:
            for line in [l.rstrip() for l in f.readlines()]:
                # search for header
                if line in self._ISA_INVESTIGATION_HEADERS:
                    if key != None:
                        tables[key] = pd.DataFrame(
                            table_lines
                        ).T  # each subtable is transposed in the i_file
                        table_lines = list()
                    key = line  # set next table key
                else:
                    tokens = line.split("\t")  # tab separated
                    table_lines.append(tokens)
        tables[key] = pd.DataFrame(
            table_lines
        ).T  # each subtable is transposed in the i_file

        # reformat each table
        def clean_quotes(string: str) -> str:
            SINGLE_OR_DOUBLE_QUOTES = "\"'"
            # don't perform on non-string elements
            if not isinstance(string, str):
                return string
            else:
                return string.lstrip(SINGLE_OR_DOUBLE_QUOTES).rstrip(
                    SINGLE_OR_DOUBLE_QUOTES
                )

        df: pd.DataFrame
        for key, df in tables.items():

            # note: as a ref, no reassign needed
            tables[key] = (
                df.rename(columns=df.iloc[0]).drop(df.index[0]).applymap(clean_quotes)
            )

        # ensure all expected subtables present
        assert set(tables.keys()) == self._ISA_INVESTIGATION_HEADERS

        self._isa_investigation_subtables = tables
        return tables


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

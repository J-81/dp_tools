""" Tools for finding specific paths 

These should rely on configuratble files to accommodate changing directory structures.
For now, the regexes will be hardcoded.

Such Locators SHOULD:
  - have find function that returns the path
  - search starting at root data directory
  - any RE patterns should be relative to the search_path
"""
import abc
import os
from pathlib import Path
from typing import List, Protocol, Tuple, runtime_checkable
import re
import logging

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)


def _rIterDir(p: Path) -> List[Path]:
    return [f for f in p.rglob("*")]


@runtime_checkable
class Locator(Protocol):
    def find(self, **kwargs) -> Path:
        """ Find a specific path """
        ...


class MultiQCDir:

    EXACT_PATH_FORMAT = os.path.join("{rel_dir}", "{mqc_label}_multiqc_report")

    def __init__(self, search_root: Path):
        self.search_root = search_root

    def find(self, rel_dir: Path, mqc_label: str) -> Path:
        found = list()
        # replace needed for regex to interpret regex escape characters AFTER interpretting python escape characters
        # (i.e. accomodate windows using the same separator as the escape char)
        pattern = self.EXACT_PATH_FORMAT.format(
            rel_dir=str(rel_dir), mqc_label=mqc_label
        ).replace("\\", r"\\")
        log.debug(
            f"Locating multiQC direcotry {type(self).__name__} with pattern: {pattern}"
        )

        putative_found = self.search_root / pattern
        assert (
            putative_found.exists()
        ), f"Target does not exist. Expected: {putative_found}. Pattern: {pattern}"
        found.append(putative_found)

        # Here is where additional validation should occur if the main found resource is a path (e.g. the RSEM stat folder)
        assert (
            len(found) == 1
        ), f"One and only find target should occur. Found: {found}. Pattern: {pattern}"
        return found[0]


class TrimmingReport:
    EXACT_PATH_FORMAT = {
        "Forward": os.path.join(
            "01-TG_Preproc",
            "Trimming_Reports",
            "{sample_name}_R1_raw.fastq.gz_trimming_report.txt",
        ),
        "Reverse": os.path.join(
            "01-TG_Preproc",
            "Trimming_Reports",
            "{sample_name}_R2_raw.fastq.gz_trimming_report.txt",
        ),
    }

    def __init__(self, search_root: Path):
        self.search_root = search_root

    def find(self, sample_name: str, type: str) -> Path:
        found = list()
        # replace needed for regex to interpret regex escape characters AFTER interpretting python escape characters
        # (i.e. accomodate windows using the same separator as the escape char)
        pattern = (
            self.EXACT_PATH_FORMAT[type]
            .format(sample_name=sample_name)
            .replace("\\", r"\\")
        )
        log.debug(f"Locating file with pattern: {pattern}")

        putative_found = self.search_root / pattern
        assert (
            putative_found.exists()
        ), f"Target does not exist. Expected: {putative_found}. Pattern: {pattern}"
        found.append(putative_found)

        # Here is where additional validation should occur if the main found resource is a path (e.g. the RSEM stat folder)
        assert (
            len(found) == 1
        ), f"One and only find target should occur. Found: {found}. Pattern: {pattern}"
        return found[0]


class FastqcReport:
    EXACT_PATH_FORMAT = {
        "Raw": {
            "Forward": os.path.join(
                "00-RawData", "FastQC_Reports", "{sample_name}_R1_raw_fastqc.{ext}"
            ),
            "Reverse": os.path.join(
                "00-RawData", "FastQC_Reports", "{sample_name}_R2_raw_fastqc.{ext}"
            ),
        },
        "Trimmed": {
            "Forward": os.path.join(
                "01-TG_Preproc",
                "FastQC_Reports",
                "{sample_name}_R1_trimmed_fastqc.{ext}",
            ),
            "Reverse": os.path.join(
                "01-TG_Preproc",
                "FastQC_Reports",
                "{sample_name}_R2_trimmed_fastqc.{ext}",
            ),
        },
    }

    def __init__(self, search_root: Path):
        self.search_root = search_root

    def find(self, sample_name: str, type: Tuple[str, str], ext: str) -> Path:
        found = list()
        # replace needed for regex to interpret regex escape characters AFTER interpretting python escape characters
        # (i.e. accomodate windows using the same separator as the escape char)
        pattern = (
            self.EXACT_PATH_FORMAT[type[0]][type[1]]
            .format(sample_name=sample_name, ext=ext)
            .replace("\\", r"\\")
        )
        log.debug(f"Locating raw fastqgz file {sample_name} with pattern: {pattern}")

        putative_found = self.search_root / pattern
        assert (
            putative_found.exists()
        ), f"Target does not exist. Expected: {putative_found}. Pattern: {pattern}"
        found.append(putative_found)

        # Here is where additional validation should occur if the main found resource is a path (e.g. the RSEM stat folder)
        assert (
            len(found) == 1
        ), f"One and only find target should occur. Found: {found}. Pattern: {pattern}"
        return found[0]


class Runsheet:
    # replace needed for regex to interpret regex escape characters AFTER interpretting python escape characters
    # (i.e. accomodate windows using the same separator as the escape char)
    EXACT_PATH_FORMAT = os.path.join(
        "Metadata",
        "AST_autogen_template_RNASeq_RCP_{datasystem_name}_RNASeq_runsheet.csv",
    ).replace("\\", r"\\")

    def __init__(self, search_root: Path):
        self.search_root = search_root

    def find(self, datasystem_name: str) -> Path:
        found = list()
        pattern = self.EXACT_PATH_FORMAT.format(datasystem_name=datasystem_name)
        log.debug(
            f"Locating {type(self).__name__} file {datasystem_name} with pattern: {pattern}"
        )

        putative_found = self.search_root / pattern
        assert (
            putative_found.exists()
        ), f"Target does not exist. Expected: {putative_found}. Pattern: {pattern}"
        found.append(putative_found)

        # Here is where additional validation should occur if the main found resource is a path (e.g. the RSEM stat folder)
        assert (
            len(found) == 1
        ), f"One and only find target should occur. Found: {found}. Pattern: {pattern}"
        return found[0]


class Fastq:
    EXACT_PATH_FORMAT = {
        "Raw": {
            "Forward": os.path.join(
                "00-RawData", "Fastq", "{sample_name}_R1_raw.fastq.gz"
            ),
            "Reverse": os.path.join(
                "00-RawData", "Fastq", "{sample_name}_R2_raw.fastq.gz"
            ),
        },
        "Trimmed": {
            "Forward": os.path.join(
                "01-TG_Preproc", "Fastq", "{sample_name}_R1_trimmed.fastq.gz"
            ),
            "Reverse": os.path.join(
                "01-TG_Preproc", "Fastq", "{sample_name}_R2_trimmed.fastq.gz"
            ),
        },
    }

    def __init__(self, search_root: Path):
        self.search_root = search_root

    def find(self, sample_name: str, type: Tuple[str, str]) -> Path:
        found = list()
        # replace needed for regex to interpret regex escape characters AFTER interpretting python escape characters
        # (i.e. accomodate windows using the same separator as the escape char)
        pattern = (
            self.EXACT_PATH_FORMAT[type[0]][type[1]]
            .format(sample_name=sample_name)
            .replace("\\", r"\\")
        )
        log.debug(f"Locating raw fastqgz file {sample_name} with pattern: {pattern}")

        putative_found = self.search_root / pattern
        assert (
            putative_found.exists()
        ), f"Target does not exist. Expected: {putative_found}. Pattern: {pattern}"
        found.append(putative_found)

        # for i, p in enumerate(_rIterDir(self.search_root)):
        #     if re.match(pattern, str(p.relative_to(self.search_root))):
        #         found.append(p)

        # Here is where additional validation should occur if the main found resource is a path (e.g. the RSEM stat folder)
        assert (
            len(found) == 1
        ), f"One and only find target should occur. Found: {found}. Pattern: {pattern}"
        return found[0]


# TODO: refactor most to this abstract class
class Exact_Path_Finder(abc.ABC):
    # replace needed for regex to interpret regex escape characters AFTER interpretting python escape characters
    # (i.e. accomodate windows using the same separator as the escape char)
    @property
    def _EXACT_PATH_FORMAT(self):
        return self.EXACT_PATH_FORMAT.replace("\\", r"\\")

    def __init__(self, search_root: Path):
        self.search_root = search_root

    def find(self, **kwargs) -> Path:
        found = list()
        pattern = self._EXACT_PATH_FORMAT.format(**kwargs)
        log.debug(f"Locating {type(self).__name__} file with pattern: {pattern}")

        putative_found = self.search_root / pattern
        assert (
            putative_found.exists()
        ), f"Target does not exist. Expected: {putative_found}. Pattern: {pattern}"
        found.append(putative_found)

        # Here is where additional validation should occur if the main found resource is a path (e.g. the RSEM stat folder)
        assert (
            len(found) == 1
        ), f"One and only find target should occur. Found: {found}. Pattern: {pattern}"
        return found[0]


class AlignedToTranscriptomeBam(Exact_Path_Finder):

    EXACT_PATH_FORMAT = os.path.join(
        "02-STAR_Alignment",
        "{sample_name}",
        "{sample_name}_Aligned.toTranscriptome.out.bam",
    )


class AlignedSortedByCoordBam(Exact_Path_Finder):

    EXACT_PATH_FORMAT = os.path.join(
        "02-STAR_Alignment",
        "{sample_name}",
        "{sample_name}_Aligned.sortedByCoord.out.bam",
    )


class AlignedSortedByCoordResortedBam(Exact_Path_Finder):

    EXACT_PATH_FORMAT = os.path.join(
        "02-STAR_Alignment",
        "{sample_name}",
        "{sample_name}_Aligned.sortedByCoord_sorted.out.bam",
    )


class AlignedSortedByCoordResortedBamIndex(Exact_Path_Finder):

    EXACT_PATH_FORMAT = os.path.join(
        "02-STAR_Alignment",
        "{sample_name}",
        "{sample_name}_Aligned.sortedByCoord_sorted.out.bam.bai",
    )


class LogFinal(Exact_Path_Finder):

    EXACT_PATH_FORMAT = os.path.join(
        "02-STAR_Alignment", "{sample_name}", "{sample_name}_Log.final.out",
    )


class LogProgress(Exact_Path_Finder):

    EXACT_PATH_FORMAT = os.path.join(
        "02-STAR_Alignment", "{sample_name}", "{sample_name}_Log.final.out",
    )


class LogFull(Exact_Path_Finder):

    EXACT_PATH_FORMAT = os.path.join(
        "02-STAR_Alignment", "{sample_name}", "{sample_name}_Log.out",
    )


class SjTab(Exact_Path_Finder):

    EXACT_PATH_FORMAT = os.path.join(
        "02-STAR_Alignment", "{sample_name}", "{sample_name}_SJ.out.tab",
    )

class GeneBodyCoverageOut(Exact_Path_Finder):

    EXACT_PATH_FORMAT = os.path.join(
        "RSeQC_Analyses", "02_geneBody_coverage", "{sample_name}"
    )

class InferExperimentOut(Exact_Path_Finder):

    EXACT_PATH_FORMAT = os.path.join(
        "RSeQC_Analyses", "03_infer_experiment", "{sample_name}_infer_expt.out"
    )

class InnerDistanceOut(Exact_Path_Finder):

    EXACT_PATH_FORMAT = os.path.join(
        "RSeQC_Analyses", "04_inner_distance", "{sample_name}"
    )

class ReadDistributionOut(Exact_Path_Finder):

    EXACT_PATH_FORMAT = os.path.join(
        "RSeQC_Analyses", "05_read_distribution", "{sample_name}_read_dist.out"
    )
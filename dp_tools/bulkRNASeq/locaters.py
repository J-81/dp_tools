""" Tools for finding specific paths 

These should rely on configuratble files to accommodate changing directory structures.
For now, the regexes will be hardcoded.

Such Locators SHOULD:
  - have find function that returns the path
  - search starting at root data directory
  - any RE patterns should be relative to the search_path
"""
from pathlib import Path
from typing import List, Protocol, Tuple, runtime_checkable
import re
import logging

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)

SEP = "/"


def _rIterDir(p: Path) -> List[Path]:
    return [f for f in p.rglob("*")]


@runtime_checkable
class Locator(Protocol):
    def find(self, **kwargs) -> Path:
        """ Find a specific path """
        ...

class MultiQCDir:
    RE_PRE_FORMAT = "{rel_dir}{SEP}{mqc_label}_multiqc_report$"

    def __init__(self, search_root: Path):
        self.search_root = search_root

    def find(self, rel_dir: Path, mqc_label: str) -> Path:
        found = list()
        pattern = self.RE_PRE_FORMAT.format(
            rel_dir=str(rel_dir), SEP=SEP, mqc_label=mqc_label
        )
        log.debug(f"Locating multiQC direcotry {type(self).__name__} with pattern: {pattern}")

        for i, p in enumerate(_rIterDir(self.search_root)):
            if re.match(pattern, str(p.relative_to(self.search_root))):
                found.append(p)

        # Here is where additional validation should occur if the main found resource is a path (e.g. the RSEM stat folder)
        assert (
            len(found) == 1
        ), f"One and only find target should occur. Found: {found}. Pattern: {pattern}"
        return found[0]

class Runsheet:
    RE_PRE_FORMAT = "Metadata{SEP}AST_autogen_template_RNASeq_RCP_{datasystem_name}_RNASeq_runsheet.csv"

    def __init__(self, search_root: Path):
        self.search_root = search_root

    def find(self, datasystem_name: str) -> Path:
        found = list()
        pattern = self.RE_PRE_FORMAT.format(
            datasystem_name=datasystem_name, SEP=SEP
        )
        log.debug(f"Locating {type(self).__name__} file {datasystem_name} with pattern: {pattern}")

        for i, p in enumerate(_rIterDir(self.search_root)):
            if re.match(pattern, str(p.relative_to(self.search_root))):
                found.append(p)

        # Here is where additional validation should occur if the main found resource is a path (e.g. the RSEM stat folder)
        assert (
            len(found) == 1
        ), f"One and only find target should occur. Found: {found}. Pattern: {pattern}"
        return found[0]

class RawFastq:
    RE_PRE_FORMAT = {
        "Raw": {
            "Forward": "00-RawData{SEP}Fastq{SEP}{sample_name}_R1_raw.fastq.gz",
            "Reverse": "00-RawData{SEP}Fastq{SEP}{sample_name}_R2_raw.fastq.gz",
        },
        "Trimmed": {
            "Forward": "01-TG_Preproc{SEP}Fastq{SEP}{sample_name}_R1_trimmed.fastq.gz",
            "Reverse": "01-TG_Preproc{SEP}Fastq{SEP}{sample_name}_R2_trimmed.fastq.gz",
        },
    }

    def __init__(self, search_root: Path):
        self.search_root = search_root

    def find(self, sample_name: str, type: Tuple[str, str]) -> Path:
        found = list()
        pattern = self.RE_PRE_FORMAT[type[0]][type[1]].format(
            sample_name=sample_name, SEP=SEP
        )
        log.debug(f"Locating raw fastqgz file {sample_name} with pattern: {pattern}")

        for i, p in enumerate(_rIterDir(self.search_root)):
            if re.match(pattern, str(p.relative_to(self.search_root))):
                found.append(p)

        # Here is where additional validation should occur if the main found resource is a path (e.g. the RSEM stat folder)
        assert (
            len(found) == 1
        ), f"One and only find target should occur. Found: {found}. Pattern: {pattern}"
        return found[0]


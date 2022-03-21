import gzip
import logging
from pathlib import Path
from typing import List
from dp_tools.components.components import RawReadsComponent
from dp_tools.core.entity_model import (
    DataDir,
    DataFile,
    TemplateDataset,
    TemplateSample,
)

log = logging.getLogger(__name__)

from dp_tools.core.check_model import Check, Flag, FlagCode


class SAMPLE_RAWREADS_0001(Check):
    id = "SAMPLE_RAWREADS_0001"
    description = (
        "Check that appropriate raw reads components exist. Also check that "
        "All datafiles associated with the components are present. "
        "For paired end studies, this means both rawForwardReads and rawReverseReads "
        "Are attached components. For single end studies, "
        "this means the rawReads component is attached. "
        "For paired end studies, confirms that forward and reverse read counts match."
    )
    flag_desc = {
        FlagCode.GREEN: "All expected raw read files present",
        FlagCode.HALT1: "Missing expected components: {missing_components}",
        FlagCode.HALT2: "Forward and reverse reads counts differ. Forward: ({forward_read_count}) Reverse: ({reverse_read_count})",
        FlagCode.DEV_HANDLED: "Searched for component, but component was not expected by entity model: {unexpected_components}",
    }

    def validate_func(self, sample: TemplateSample) -> Flag:
        # assume passing unless a flag condition arises
        code = FlagCode.GREEN

        # set branching informative parameters based on layout
        if sample.dataset.metadata.paired_end:
            expected_components = ["rawForwardReads", "rawReverseReads"]
            check_read_parity = True
        else:
            expected_components = ["rawReads"]
            check_read_parity = False

        missing_components = list()
        unexpected_components = list()
        for expected_component in expected_components:
            component = getattr(sample, expected_component, None)
            if component == None:
                unexpected_components.append(expected_component)
            if not isinstance(component, RawReadsComponent):
                missing_components.append(expected_component)

        if unexpected_components:
            code = FlagCode.DEV_HANDLED

        if missing_components:
            code = FlagCode.HALT1

        # check parity
        if all([check_read_parity, code == FlagCode.GREEN]):
            if (
                not sample.rawForwardReads.mqcData["FastQC"]["General_Stats"][
                    "Total Sequences"
                ]
                == sample.rawReverseReads.mqcData["FastQC"]["General_Stats"][
                    "Total Sequences"
                ]
            ):
                code = FlagCode.HALT2

        return Flag(
            check=self,
            code=code,
            message_args={
                "missing_components": missing_components,
                "forward_read_count": sample.rawForwardReads.mqcData["FastQC"][
                    "General_Stats"
                ]["Total Sequences"]
                if code == FlagCode.HALT2
                else None,
                "reverse_read_count": sample.rawReverseReads.mqcData["FastQC"][
                    "General_Stats"
                ]["Total Sequences"]
                if code == FlagCode.HALT2
                else None,
            },
        )

class SAMPLE_TRIMREADS_0001(SAMPLE_RAWREADS_0001):
    ...


class COMPONENT_RAWREADS_0001(Check):
    id = "COMPONENT_RAWREADS_0001"
    config = {
        "lines_to_check": 200_000_000,
        # attributes names
        "expected_data_files": [
            "fastqGZ",
            "multiQCDir",
            "fastqcReportHTML",
            "fastqcReportZIP",
        ],
    }
    description = (
        "Confirms that all read components (e.g. rawForwardReads, trimmedReads) should include the following: "
        "Datafiles of the format: {expected_data_files} related to the reads component. "
        "Additionally, the following checks are performed for each file type: \n"
        "\tfastq.gz: First {lines_to_check} lines are checked for correct format. "
    )
    flag_desc = {
        FlagCode.GREEN: "Component passes all validation requirements.",
        FlagCode.HALT1: "Missing expected files: {missing_files}",
        FlagCode.HALT2: "Fastq.gz file has issues on lines: {lines_with_issues}",
        FlagCode.HALT3: "Corrupted Fastq.gz file suspected, last line number encountered: {last_line_checked}",
    }

    def validate_func(self: Check, component) -> Flag:
        """ Checks fastq lines for expected header content
        Note: Example of header from GLDS-194
        |  ``@J00113:376:HMJMYBBXX:3:1101:26666:1244 1:N:0:NCGCTCGA\n``
        This also assumes the fastq file does NOT split sequence or quality lines
        for any read
        :param component: A ReadsComponent
        """
        # assume passing first
        # overwrite if flag conditions met
        code = FlagCode.GREEN

        # Subcheck: 1 ( can trigger HALT1 )
        # check if expected files exist first
        missing_files: List[Path] = list()
        lines_with_issues: List[int] = list()
        i = 0

        for expected_file in self.config["expected_data_files"]:
            try:
                # check the attribute is exists and is of the proper type
                assert any(
                    [
                        isinstance(getattr(component, expected_file), DataFile),
                        isinstance(getattr(component, expected_file), DataDir),
                    ]
                )
                # check the path exists
                assert getattr(component, expected_file).path.exists()
            except AssertionError:
                code = FlagCode.HALT1
                missing_files.append(expected_file)

        # check if exiting makes sense before next checks
        if code != FlagCode.GREEN:
            return Flag(
                check=self,
                code=code,
                message_args={
                    "lines_with_issues": lines_with_issues,
                    "last_line_checked": i,
                    "missing_files": missing_files,
                },
            )

        # subcheck: 2 ( can trigger HALT2,HALT3 )
        # check fastq.gz file looks correct
        file = component.fastqGZ.path
        count_lines_to_check = self.config["lines_to_check"]

        if count_lines_to_check == -1:
            count_lines_to_check = float("inf")

        # truncated files raise EOFError
        # catch this as HALT3
        try:
            with gzip.open(file, "rb") as f:
                for i, line in enumerate(f):
                    # checks if lines counted equals the limit input
                    if i + 1 == count_lines_to_check:
                        log.debug(
                            f"Reached {count_lines_to_check} lines, ending line check"
                        )
                        break

                    line = line.decode()
                    # every fourth line should be an identifier
                    expected_identifier_line = i % 4 == 0
                    # check if line is actually an identifier line
                    if expected_identifier_line and line[0] != "@":
                        lines_with_issues.append(i + 1)
                    # update every 20,000,000 reads
                    if i % 20_000_000 == 0:
                        log.debug(f"Checked {i} lines for {file}")
                        pass
            if not len(lines_with_issues) == 0:
                code = FlagCode.HALT2
        except (EOFError, gzip.BadGzipFile):
            code = FlagCode.HALT3

        # return flag
        return Flag(
            check=self,
            code=code,
            message_args={
                "lines_with_issues": lines_with_issues,
                "last_line_checked": i,
                "missing_files": missing_files,
            },
        )


class COMPONENT_TRIMREADS_0001(Check):
    id = "COMPONENT_TRIMREADS_0001"
    config = {
        "lines_to_check": 200_000_000,
        "expected_data_files": [
            "fastqGZ",
            "multiQCDir",
            "fastqcReportHTML",
            "fastqcReportZIP",
            "trimmingReportTXT",
        ],
    }


class DATASET_RAWREADS_0001(Check):
    id = "DATASET_RAWREADS_0001"
    config = {
        "lines_to_check": 200_000_000,
        # attributes names
        "expected_data_files": [
            "fastqGZ",
            "multiQCDir",
            "fastqcReportHTML",
            "fastqcReportZIP",
        ],
    }
    description = (
        "Performs a validation across all samples' raw reads "
        "Specifically, checks for any sample-wise outliers in the following metrics: "
        "\n\t fastqGZ file size, read counts. "
        "\Also, checks for any sample-wise outliers in the following metrics with the read sets: "
        "\n\t read counts. "
        "Additionally, the following checks are performed for each file type: \n"
        "\tfastq.gz: First {lines_to_check} lines are checked for correct format. "
    )
    flag_desc = {
        FlagCode.GREEN: "Component passes all validation requirements.",
        FlagCode.HALT1: "Missing expected files: {missing_files}",
        FlagCode.HALT2: "Fastq.gz file has issues on lines: {lines_with_issues}",
        FlagCode.HALT3: "Corrupted Fastq.gz file suspected, last line number encountered: {last_line_checked}",
    }

    def validate(self: Check, dataset: TemplateDataset) -> Flag:
        return Flag(code=FlagCode.GREEN, check=self)


class DATASET_TRIMREADS_0001(DATASET_RAWREADS_0001):
    id = "DATASET_TRIMREADS_0001"
    config = {
        "lines_to_check": 200_000_000,
        "expected_data_files": [
            "fastqGZ",
            "multiQCDir",
            "fastqcReportHTML",
            "fastqcReportZIP",
            "trimmingReportTXT",
        ],
    }

import gzip
import logging
from dp_tools.components.components import RawReadsComponent

log = logging.getLogger(__name__)

from dp_tools.core.check_model import Check, Flag, FlagCode


def _validate_func_SAMPLE_RAWREADS_0001(self: Check, sample) -> Flag:
    if sample.dataset.metadata.paired_end:
        expected_components = ["rawForwardReads", "rawReverseReads"]
    else:
        expected_components = ["rawReads"]

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
    else:
        code = FlagCode.GREEN
    return Flag(
        check=self, code=code, message_args={"missing_components": missing_components}
    )


SAMPLE_RAWREADS_0001 = Check(
    id="SAMPLE_RAWREADS_0001",
    description=(
        "Check that appropriate raw reads components exist. Also check that "
        "All datafiles associated with the components are present. "
        "For paired end studies, this means both rawForwardReads and rawReverseReads "
        "Are attached components. For single end studies, "
        "this means the rawReads component is attached. "
    ),
    flag_desc={
        FlagCode.GREEN: "All expected raw read files present",
        FlagCode.HALT1: "Missing expected components: {missing_components}",
        FlagCode.DEV_HANDLED: "Searched for component, but component was not expected by entity model: {unexpected_components}",
    },
    validate_func=_validate_func_SAMPLE_RAWREADS_0001,
)


def _validate_func_COMPONENT_READS(self: Check, component) -> Flag:
    """ Checks fastq lines for expected header content
    Note: Example of header from GLDS-194
    |  ``@J00113:376:HMJMYBBXX:3:1101:26666:1244 1:N:0:NCGCTCGA\n``
    This also assumes the fastq file does NOT split sequence or quality lines
    for any read
    :param component: A ReadsComponent
    """
    file = component
    count_lines_to_check = self.config["lines_to_check"]
    expected_files = self.config["expected_data_files"]

    if count_lines_to_check == -1:
        count_lines_to_check = float("inf")

    lines_with_issues = list()

    # check for expected files

    # truncated files raise EOFError
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
        if len(lines_with_issues) != 0:
            code = FlagCode.GREEN
        else:
            code = FlagCode.HALT2
    except EOFError:
        code = FlagCode.HALT3
    return Flag(
        check=self, code=code, message_args={"lines_with_issues": lines_with_issues}
    )


COMPONENT_RAWREADS_0001 = Check(
    config={
        "lines_to_check": 200_000_000,
        "expected_data_files": [
            "fastqGZ",
            "multiQCDir",
            "fastqcReportHTML",
            "fastqcReportZIP",
        ],
    },
    id="COMPONENT_RAWREADS_0001",
    description=(
        "Confirms that all read components (e.g. rawForwardReads, trimmedReads) should include the following: "
        "Datafiles of the format: {expected_data_files} related to the reads component. "
        "Additionally, the following checks are performed for each file type: \n"
        "\tfastq.gz: First {lines_to_check} lines are checked for correct format. "
    ),
    flag_desc={
        FlagCode.GREEN: "Component passes all validation requirements.",
        FlagCode.HALT1: "Missing expected files: {missing_files}",
        FlagCode.HALT2: "Fastq.gz file has issues on lines: {lines_with_issues}",
        FlagCode.HALT3: "Corrupted Fastq.gz file suspected: {last_line_checked}",
    },
    validate_func=_validate_func_COMPONENT_READS,
)

COMPONENT_TRIMREADS_0001 = COMPONENT_RAWREADS_0001.copy_with_new_config(
    id="COMPONENT_TRIMREADS_0001",
    config={
        "lines_to_check": 200_000_000,
        "expected_data_files": [
            "fastqGZ",
            "multiQCDir",
            "fastqcReportHTML",
            "fastqcReportZIP",
            "trimmingReportTXT",
        ],
    },
)

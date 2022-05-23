from collections import defaultdict
import copy
import gzip
import logging
import math
from pathlib import Path
import subprocess
from typing import Dict, Union
from dp_tools.bulkRNASeq.entity import BulkRNASeqDataset

import pandas as pd
from dp_tools.components.components import RawReadsComponent
from dp_tools.core.entity_model import (
    BaseComponent,
    ModuleLevelMQC,
    TemplateComponent,
    TemplateDataset,
)

log = logging.getLogger(__name__)

from dp_tools.core.check_model import FlagCode, FlagEntry, FlagEntryWithOutliers

# adapted from reference: https://stackoverflow.com/questions/56048627/round-floats-in-a-nested-dictionary-recursively
# used to round values for easier to read messages
def formatfloat(x):
    return "%.3g" % float(x)


def pformat(original_dictionary, function):
    dictionary = copy.deepcopy(
        original_dictionary
    )  # we don't want to override original values
    if isinstance(dictionary, dict):
        new_dict = dict()
        for k, v in dictionary.items():
            new_dict[k] = function(v) if isinstance(v, float) else pformat(v, function)
        return new_dict
    return dictionary


def convert_nan_to_zero(input: Dict[str, Union[float, int]]) -> Dict:
    """Convert any Nan into zero"""
    output = dict()
    for key, value in input.items():
        output[key] = value if not math.isnan(value) else 0
    return output


## Functions that use the following syntax to merge values from general stats:
# "stat1 + stat2" should search and sum the stats
def stat_string_to_value(stat_string: str, mqcData: ModuleLevelMQC) -> float:
    """ "stat1 + stat2" should search and sum the stats"""
    sum = float(0)
    direct_keys = stat_string.split(" + ")
    for direct_key in direct_keys:
        print(direct_key)
        sum += mqcData["General_Stats"][direct_key]
    return sum


## Dataframe and Series specific helper functions
def nonNull(df: pd.DataFrame) -> bool:
    # negation since it checks if any are null
    return ~df.isnull().any(axis=None)


def nonNegative(df: pd.DataFrame) -> bool:
    """This ignores null values, use nonNull to validate that condition"""
    return ((df >= 0) | (df.isnull())).all(axis=None)


def onlyAllowedValues(df: pd.DataFrame, allowed_values: list) -> bool:
    """This ignores null values, use nonNull to validate that condition"""
    return ((df.isin(allowed_values)) | (df.isnull())).all(axis=None)


def check_forward_and_reverse_reads_counts_match(
    fwd_reads: RawReadsComponent, rev_reads: RawReadsComponent
) -> FlagEntry:
    # data specific preprocess
    count_fwd_reads = fwd_reads.mqcData["FastQC"]["General_Stats"]["total_sequences"]
    count_rev_reads = rev_reads.mqcData["FastQC"]["General_Stats"]["total_sequences"]

    # check logic
    if count_fwd_reads == count_rev_reads:
        code = FlagCode.GREEN
        message = (
            f"Forward and reverse read counts match at "
            f"{int(count_rev_reads)} sequences "
        )
    else:
        code = FlagCode.HALT
        message = (
            f"Forward and reverse read counts do not "
            f"match: forward_Count:{int(count_fwd_reads)}, "
            f"reverse_Count:{int(count_rev_reads)}"
        )

    return {"code": code, "message": message}


def check_file_exists(file: Path) -> FlagEntry:
    # check logic
    if file.is_file():
        code = FlagCode.GREEN
        message = f"File exists: {file.name} "
    else:
        code = FlagCode.HALT
        message = f"Missing file: {file.name} expected at {str(file)} "

    return {"code": code, "message": message}


def check_fastqgz_file_contents(file: Path, count_lines_to_check: int) -> FlagEntry:
    """Check fastqgz by:
    1. Decompressing as a stream of lines.
    2. Affirming expected headers (every 4th line) look correct.

    :param file: Input fastqGZ file path
    :type file: Path
    :param count_lines_to_check: Maximum number of lines to check. Setting this to a negative value will remove the limit
    :type count_lines_to_check: int
    :return: A required fields-only flag entry dictionary
    :rtype: FlagEntry
    """

    lines_with_issues: list[int] = list()

    # check logic
    # truncated files raise EOFError
    # catch this as HALT3
    try:
        with gzip.open(file, "rb") as f:
            for i, byte_line in enumerate(f):
                # checks if lines counted equals the limit input
                if i + 1 == count_lines_to_check:
                    log.debug(
                        f"Reached {count_lines_to_check} lines, ending line check"
                    )
                    break

                line = byte_line.decode()
                # every fourth line should be an identifier
                expected_identifier_line = i % 4 == 0
                # check if line is actually an identifier line
                if expected_identifier_line and line[0] != "@":
                    lines_with_issues.append(i + 1)
                # update every 2,000,000 reads
                if i % 2_000_000 == 0:
                    log.debug(f"Checked {i} lines for {file}")
                    pass

        if not len(lines_with_issues) == 0:
            code = FlagCode.HALT
            message = (
                f"Following decompressed fastqGZ lines have issues: {lines_with_issues}"
            )
        else:
            code = FlagCode.GREEN
            message = f"First {count_lines_to_check} lines checked found no issues."
    except (EOFError, gzip.BadGzipFile):
        code = FlagCode.HALT
        message = (
            f"Error during decompression, likely a compression or truncation issue."
        )

    return {"code": code, "message": message}


def check_bam_file_integrity(file: Path, samtools_bin: Path) -> FlagEntry:
    """Uses http://www.htslib.org/doc/samtools-quickcheck.html"""
    # data specific preprocess

    # check logic
    output = subprocess.run(
        [str(samtools_bin), "quickcheck", "-v", str(file)], capture_output=True
    )
    stdout_string = output.stdout.decode()
    if stdout_string == "":
        code = FlagCode.GREEN
        message = f"Samtools quickcheck raised no issues"
    else:
        code = FlagCode.HALT
        message = (
            f"Samtools quickcheck failed on this file with output: {stdout_string}"
        )
    return {"code": code, "message": message}


def check_thresholds(
    component: TemplateComponent, mqc_key: str, stat_string: str, thresholds: list[dict]
) -> FlagEntry:
    # data specific preprocess
    value = stat_string_to_value(stat_string, component.mqcData[mqc_key])

    # check logic
    # Assuming GREEN unless reassigned
    code = FlagCode.GREEN
    for threshold in thresholds:
        match threshold["type"]:
            case "lower":
                if value < threshold["value"]:
                    code = (
                        FlagCode[threshold["code"]]
                        if code < FlagCode[threshold["code"]]
                        else code
                    )

    if code == FlagCode.GREEN:
        message = f"Value: ({value}) did not breech any configured thresholds"
    else:
        message = f"Value: ({value}) breeched configured thresholds"
    return {"code": code, "message": message}


def check_metadata_attributes_exist(
    dataset: BulkRNASeqDataset, expected_attrs: list[str]
) -> FlagEntry:
    # data specific preprocess
    # set up tracker for expected attributes values
    tracked_metadata = dict()
    # and a tracker for missing attributes
    missing_metadata_fields = list()

    for attr in expected_attrs:
        attr_value = getattr(dataset.metadata, attr, None)
        if attr_value != None:
            tracked_metadata[attr] = attr_value
        else:
            missing_metadata_fields.append(attr)

    # check if any missing_metadata_fields are present
    # check logic
    if not missing_metadata_fields:
        code = FlagCode.GREEN
        message = f"All expected metadata values found: {tracked_metadata}"
    else:
        code = FlagCode.HALT
        message = f"Missing dataset metadata (source from Runsheet): {missing_metadata_fields}"
    return {"code": code, "message": message}


def check_for_outliers(
    dataset: TemplateDataset,
    sample_component: str,
    mqc_module: str,
    mqc_plot: str,
    mqc_keys: list[str],
    thresholds: list[dict],
) -> FlagEntryWithOutliers:
    # assume code is GREEN until outliers detected
    code = FlagCode.GREEN
    # dataframe extraction
    df: pd.DataFrame = dataset.getMQCDataFrame(
        sample_component=sample_component, mqc_module=mqc_module, mqc_plot=mqc_plot
    )

    def default_to_regular(d):
        if isinstance(d, defaultdict):
            d = {k: default_to_regular(v) for k, v in d.items()}
        return d

    # track for outliers
    outliers: dict[str, dict[str, dict[str, str]]] = defaultdict(
        lambda: defaultdict(dict)
    )

    # override if mqc_keys is a special value
    if mqc_keys == ["_ALL"]:
        mqc_keys = df.columns

    for mqc_key in mqc_keys:
        for threshold in thresholds:
            if threshold["middle_fcn"] == "mean":
                middle = df[mqc_key].mean()
            elif threshold["middle_fcn"] == "median":
                middle = df[mqc_key].median()
            else:
                raise ValueError(
                    f"Cannot compute middle from supplied middle_fcn name: {threshold['middle_fcn']}. Must supply either 'median' or 'mean'"
                )

            # bail if standard deviation == 0
            # e.g. if all values are identical (and thus has no outliers)
            if df[mqc_key].std() == 0:
                continue

            # compute difference
            df_diffs = df[mqc_key] - middle

            # compute as number of standard deviations
            df_diffs_in_std = df_diffs / df[mqc_key].std()

            # add to outlier tracker if over the threshold
            for key, value in df_diffs_in_std.iteritems():
                # if an outlier
                if abs(value) > threshold["stdev_threshold"]:
                    # track it
                    outliers[key][mqc_module][mqc_key] = value
                    # elevate code if current code is lower severity
                    if code < FlagCode[threshold["code"]]:
                        code = FlagCode[threshold["code"]]
    # check logic
    if code == FlagCode.GREEN:
        message = f"No outliers found for {mqc_keys} in {mqc_plot} part of {mqc_module} multiQC module"
    else:
        message = (
            f"Outliers found in {mqc_module} multiQC module as follows: {outliers}"
        )
    return {"code": code, "message": message, "outliers": default_to_regular(outliers)}


def _check_expected_files_exist(
    input_dir: Path, expected_extensions: list[str], parent_dir_is_filename: bool = True
):
    if parent_dir_is_filename:
        fname = input_dir.name
    expected_files = [input_dir / f"{fname}{ext}" for ext in expected_extensions]
    missing_files = list()
    for expected_file in expected_files:
        if not expected_file.is_file():
            missing_files.append(str(expected_file))

    expected_file_str = [str(f) for f in expected_files]
    return missing_files, expected_file_str


def check_genebody_coverage_output(input_dir: Path):
    EXPECTED_EXTENSIONS = [
        ".geneBodyCoverage.r",
        ".geneBodyCoverage.txt",
        ".geneBodyCoverage.curves.pdf",
    ]

    missing_files, expected_file_str = _check_expected_files_exist(
        input_dir, expected_extensions=EXPECTED_EXTENSIONS
    )

    if not missing_files:
        code = FlagCode.GREEN
        message = f"All output from geneBody coverage found: {expected_file_str}"
    else:
        code = FlagCode.HALT
        message = f"Missing output from geneBody coverage: {missing_files}. Expected: {expected_file_str}"
    return {"code": code, "message": message}


def check_inner_distance_output(input_dir: Path):
    EXPECTED_EXTENSIONS = [
        ".inner_distance_plot.r",
        ".inner_distance_freq.txt",
        ".inner_distance.txt",
        ".inner_distance_plot.pdf",
    ]

    missing_files, expected_file_str = _check_expected_files_exist(
        input_dir, expected_extensions=EXPECTED_EXTENSIONS
    )

    if not missing_files:
        code = FlagCode.GREEN
        message = f"All output from inner distance found: {expected_file_str}"
    else:
        code = FlagCode.HALT
        message = f"Missing output from inner distance: {missing_files}. Expected: {expected_file_str}"
    return {"code": code, "message": message}


def check_strandedness_assessable_from_infer_experiment(
    dataset: BulkRNASeqDataset,
    stranded_assessment_range: dict[str, float],
    unstranded_assessment_range: dict[str, float],
    valid_dominant_strandedness_assessments: list[str],
) -> FlagEntry:
    # data specific preprocess
    def get_median_strandedness(
        dataset: TemplateDataset,
    ) -> dict[str, float]:
        df = dataset.getMQCDataFrame(
            sample_component="rSeQCAnalysis",
            mqc_module="RSeQC",
            mqc_plot="Infer experiment",
        ).fillna(
            0
        )  # Nan is a zero for this MultiQC table

        median_strandedness = df.median().to_dict()

        return median_strandedness

    median_strandedness = get_median_strandedness(dataset)

    # check if dominant assessment is valid
    strand_assessment: str = max(
        median_strandedness, key=lambda k: median_strandedness[k]
    )

    # flag based on thresholds
    assessment_value: float = median_strandedness[strand_assessment]

    is_stranded: bool = (
        stranded_assessment_range["max"]
        > assessment_value
        > stranded_assessment_range["min"]
    )
    is_unstranded: bool = (
        unstranded_assessment_range["max"]
        > assessment_value
        > unstranded_assessment_range["min"]
    )

    def determine_samples_outside_range(
        dataset: TemplateDataset, min: float, max: float
    ) -> list[str]:
        df = dataset.getMQCDataFrame(
            sample_component="rSeQCAnalysis",
            mqc_module="RSeQC",
            mqc_plot="Infer experiment",
        ).fillna(
            0
        )  # Nan is a zero for this MultiQC table

        return df.index[df[strand_assessment].between(min, max) == False].to_list()

    # Catalog and flag any samples outside of range
    # flags based on samples that are out of the assessment range
    samples_outside_range: list[str]
    if is_stranded:
        samples_outside_range = determine_samples_outside_range(
            dataset,
            stranded_assessment_range["min"],
            stranded_assessment_range["max"],
        )
    elif is_unstranded:
        samples_outside_range = determine_samples_outside_range(
            dataset,
            unstranded_assessment_range["min"],
            unstranded_assessment_range["max"],
        )
    else:  # this means that the strandedness is ambiguous
        samples_outside_range = list()

    # check logic
    if strand_assessment not in valid_dominant_strandedness_assessments:
        code = FlagCode.HALT
        message = f"Dominant strandedness [{strand_assessment} (median:{assessment_value:.2f})] is invalid for processing. Valid assessments: {valid_dominant_strandedness_assessments}"
    elif not samples_outside_range and any([is_stranded, is_unstranded]):
        code = FlagCode.GREEN
        message = f"Dominant strandedness [{strand_assessment} (median:{assessment_value:.2f})] assessed with no individual samples outside the assessment range"
    elif samples_outside_range and any([is_stranded, is_unstranded]):
        code = FlagCode.RED
        message = f"Dominant strandedness [{strand_assessment} (median:{assessment_value:.2f})] assessed with samples outside the assessment range: {samples_outside_range}"
    else:
        code = FlagCode.HALT
        message = (
            f"Dominant strandedness [{strand_assessment} (median:{assessment_value:.2f})] is ambiguous due to being inside range "
            f"({stranded_assessment_range['min']}-{unstranded_assessment_range['max']})"
        )

    return {"code": code, "message": message}


def check_rsem_counts_and_unnormalized_tables_parity(
    rsem_table_path: Path, deseq2_table_path: Path
) -> FlagEntry:
    # data specific preprocess
    df_rsem = pd.read_csv(rsem_table_path)
    df_deseq2 = pd.read_csv(deseq2_table_path)

    # return halt flag if column labels not conserved
    if not set(df_deseq2.columns) == set(df_rsem.columns):
        unique_to_deseq2 = set(df_deseq2.columns) - set(df_rsem.columns)
        unique_to_rsem = set(df_rsem.columns) - set(df_deseq2.columns)
        return {
            "code": FlagCode.HALT,
            "message": f"Columns do not match: unique to rsem: {unique_to_rsem}. unique to deseq2: {unique_to_deseq2}.",
        }

    # rearrange columns to the same order
    df_deseq2 = df_deseq2[df_rsem.columns]

    # check logic
    if df_deseq2.equals(df_rsem):
        code = FlagCode.GREEN
        message = f"Tables of unnormalized counts match."
    else:
        code = FlagCode.HALT
        message = (
            f"Tables of unnormalized counts have same columns but values do not match."
        )
    return {"code": code, "message": message}

from collections import defaultdict
import copy
import enum
import gzip
import itertools
import logging
import math
from pathlib import Path
from statistics import mean
import string
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


def r_style_make_names(s: str) -> str:
    """Recreates R's make.names function for individual strings.
    This function is often used to create syntactically valid names in R which are then saved in R outputs.
    Source: https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/make.names

    Args:
        s (str): A string to convert

    Returns:
        str: A string converted in the same way as R's make.names function
    """
    VALID_CHARACTERS = string.ascii_letters + string.digits + "."
    REPLACEMENT_CHAR = "."
    new_string_chars = list()
    for char in s:
        if char in VALID_CHARACTERS:
            new_string_chars.append(char)
        else:
            new_string_chars.append(REPLACEMENT_CHAR)
    return "".join(new_string_chars)


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


def check_aggregate_star_unnormalized_counts_table_values_against_samplewise_tables(
    unnormalizedCountTable: Path, samplewise_tables: dict[str, Path]
) -> FlagEntry:
    STAR_COUNT_MODES = ["unstranded", "sense", "antisense"]
    # data specific preprocess
    df_agg = pd.read_csv(unnormalizedCountTable, index_col=0)

    # based on which column matches the first entry
    # all columns must match with the same strand column
    strand_assessment: str = None  # type: ignore
    samples_with_issues: dict[str, list[str]] = {
        "Not in aggregate table": list(),
        "Sample counts mismatch": list(),
    }
    for sample, path in samplewise_tables.items():
        # check if samples exist as a column
        if sample not in df_agg:
            samples_with_issues["Not in aggregate table"].append(sample)
            break

        # load
        df_samp = pd.read_csv(
            path, sep="\t", names=STAR_COUNT_MODES, index_col=0
        ).filter(
            regex="^(?!N_.*).*", axis="rows"
        )  # filter out N_* entries

        # check if the values match for any of the count modes
        #   unstranded, sense, antisense
        # for remaining samples, only check the match for the first count mode
        for count_mode in STAR_COUNT_MODES:
            # make sure to sort indicies
            if df_agg[sample].sort_index().equals(df_samp[count_mode].sort_index()):
                # assign strand assessment if first sample
                if strand_assessment is None:
                    strand_assessment = count_mode

                if strand_assessment == count_mode:
                    # no issues found (i.e. counts match with a consistent count mode column), break out
                    break
        else:  # no break
            samples_with_issues["Sample counts mismatch"].append(sample)

    # check logic
    if not any([issue_type for issue_type in samples_with_issues.values()]):
        code = FlagCode.GREEN
        message = (
            f"All samples accounted for and with matching counts "
            f"between samplewise and aggregate table using strand assessment: '{strand_assessment}'"
        )
    else:
        code = FlagCode.HALT
        message = f"Identified issues: {samples_with_issues}"
    return {"code": code, "message": message}


def check_aggregate_rsem_unnormalized_counts_table_values_against_samplewise_tables(
    unnormalizedCountTable: Path, samplewise_tables: dict[str, Path]
) -> FlagEntry:
    # data specific preprocess
    df_agg = pd.read_csv(unnormalizedCountTable, index_col=0)

    # based on which column matches the first entry
    samples_with_issues: dict[str, list[str]] = {
        "Not in aggregate table": list(),
        "Sample counts mismatch": list(),
    }
    for sample, path in samplewise_tables.items():
        # check if samples exist as a column
        if sample not in df_agg:
            samples_with_issues["Not in aggregate table"].append(sample)
            break

        # load
        df_samp = pd.read_csv(path, sep="\t", index_col=0)  # filter out N_* entries

        # check if values match
        if (
            not df_agg[sample]
            .sort_index()
            .equals(df_samp["expected_count"].sort_index())
        ):
            samples_with_issues["Sample counts mismatch"].append(sample)

    # check logic
    if not any([issue_type for issue_type in samples_with_issues.values()]):
        code = FlagCode.GREEN
        message = f"All samples accounted for and with matching counts between samplewise and aggregate table"
    else:
        code = FlagCode.HALT
        message = f"Identified issues: {samples_with_issues}"
    return {"code": code, "message": message}


def check_sample_table_for_all_samples(runsheet: Path, sampleTable: Path) -> FlagEntry:
    """Check the sample table includes all samples as denoted in the runsheet.

    Args:
        runsheet (Path): csv file used for processing, the index denotes all samples
        sampleTable (Path): csv file that pairs each sample with resolved experimental group (called condition within the table)

    Returns:
        FlagEntry: A check result
    """
    # data specific preprocess
    df_rs = pd.read_csv(runsheet, index_col="Sample Name").sort_index()
    df_sample = pd.read_csv(sampleTable, index_col=0).sort_index()

    extra_samples: dict[str, set[str]] = {
        "unique_to_runsheet": set(df_rs.index) - set(df_sample.index),
        "unique_to_sampleTable": set(df_sample.index) - set(df_rs.index),
    }

    # check logic
    if not any([entry for entry in extra_samples.values()]):
        code = FlagCode.GREEN
        message = f"All samples matching"
    else:
        code = FlagCode.HALT
        message = f"Samples mismatched: {[f'{entry}:{v}' for entry, v  in extra_samples.items() if v]}"
    return {"code": code, "message": message}


class GroupFormatting(enum.Enum):
    r_make_names = enum.auto()
    ampersand_join = enum.auto()


def utils_runsheet_to_expected_groups(
    runsheet: Path,
    formatting: GroupFormatting = GroupFormatting.ampersand_join,
    map_to_lists: bool = False,
) -> dict[str, str]:
    df_rs = (
        pd.read_csv(runsheet, index_col="Sample Name")
        .filter(regex="^Factor Value\[.*\]")
        .sort_index()
    )  # using only Factor Value columns

    match formatting:
        case GroupFormatting.r_make_names:
            expected_conditions_based_on_runsheet = (
                df_rs.apply(lambda x: "...".join(x), axis="columns")
                .apply(r_style_make_names)  # join factors with '...'
                .to_dict()
            )  # reformat entire group in the R style
        case GroupFormatting.ampersand_join:
            expected_conditions_based_on_runsheet = df_rs.apply(
                lambda x: f"({' & '.join(x)})", axis="columns"
            ).to_dict()
        case _:
            raise ValueError(
                f"Formatting method invalid, must be one of the following: {list(GroupFormatting)}"
            )

    # convert from {sample: group} dict
    #         to {group: [samples]} dict
    if map_to_lists:
        unique_groups = set(expected_conditions_based_on_runsheet.values())
        reformatted_dict: dict[str, list[str]] = dict()
        for query_group in unique_groups:
            reformatted_dict[query_group] = [
                sample
                for sample, group in expected_conditions_based_on_runsheet.items()
                if group == query_group
            ]
        expected_conditions_based_on_runsheet = reformatted_dict

    return expected_conditions_based_on_runsheet


def check_sample_table_for_correct_group_assignments(
    runsheet: Path, sampleTable: Path
) -> FlagEntry:
    """Check the sample table is assigned to the correct experimental group.
    An experimental group is defined by the Factor Value columns found in the runsheet.

    Args:
        runsheet (Path): csv file used for processing, includes metadata used for experimental group designation
        sampleTable (Path): csv file that pairs each sample with resolved experimental group (called condition within the table)

    Returns:
        FlagEntry: A check result
    """
    # data specific preprocess
    df_rs = (
        pd.read_csv(runsheet, index_col="Sample Name")
        .filter(regex="^Factor Value\[.*\]")
        .sort_index()
    )  # using only Factor Value columns
    df_sample = pd.read_csv(sampleTable, index_col=0).sort_index()

    # TODO: refactor with utils_runsheet_to_expected_groups
    expected_conditions_based_on_runsheet = df_rs.apply(
        lambda x: "...".join(x), axis="columns"
    ).apply(  # join factors with '...'
        r_style_make_names
    )  # reformat entire group in the R style

    mismatched_rows = expected_conditions_based_on_runsheet != df_sample["condition"]

    # check logic
    if not any(mismatched_rows):
        code = FlagCode.GREEN
        message = f"Conditions are formatted and assigned correctly based on runsheet"
    else:
        code = FlagCode.HALT
        mismatch_description = (
            df_sample[mismatched_rows]["condition"]
            + " <--SAMPLETABLE : RUNSHEET--> "
            + expected_conditions_based_on_runsheet[mismatched_rows]
        ).to_dict()
        message = f"Mismatch in expected conditions based on runsheet for these rows: {mismatch_description}"
    return {"code": code, "message": message}


def check_contrasts_table_headers(contrasts_table: Path, runsheet: Path) -> FlagEntry:
    # data specific preprocess
    expected_groups = utils_runsheet_to_expected_groups(runsheet, map_to_lists=True)
    expected_comparisons = [
        "v".join(paired_groups)
        for paired_groups in itertools.permutations(expected_groups, 2)
    ]
    df_contrasts = pd.read_csv(contrasts_table, index_col=0)

    # check logic
    differences = set(expected_comparisons).symmetric_difference(
        set(df_contrasts.columns)
    )
    if not differences:
        code = FlagCode.GREEN
        message = f"Contrasts header includes expected comparisons as determined runsheet Factor Value Columns: {set(expected_comparisons)}"
    else:
        code = FlagCode.HALT
        message = f"Contrasts header does not match expected comparisons as determined runsheet Factor Value Columns: {differences}"
    return {"code": code, "message": message}


def check_contrasts_table_rows(contrasts_table: Path, **_) -> FlagEntry:
    # data specific preprocess
    df_contrasts = pd.read_csv(contrasts_table, index_col=0)

    def _get_groups_from_comparisions(s: str) -> set[str]:
        """Converts '(G1)v(G2)'
        into G1...G2 where G1 and G2 are renamed as per the r make names function

        Args:
            s (str): Input that fits this format: '(G1)v(G2)'

        Returns:
            str: Reformatted string
        """
        g1, g2 = s.split("v")
        # remove parens and reformat with r make names style
        g1 = r_style_make_names(g1[1:-1].replace(" & ", "..."))
        g2 = r_style_make_names(g2[1:-1].replace(" & ", "..."))
        return {g1, g2}

    bad_columns: list[str] = list()
    for (col_name, col_series) in df_contrasts.iteritems():
        expected_values = _get_groups_from_comparisions(col_name)
        if not expected_values == set(col_series):
            bad_columns.append(col_name)

    # check logic
    if not bad_columns:
        code = FlagCode.GREEN
        message = f"Contrasts column and rows match expected formatting"
    else:
        code = FlagCode.HALT
        message = f"Contrasts columns {bad_columns} have unexpected values"
    return {"code": code, "message": message}


def check_dge_table_annotation_columns_exist(
    dge_table: Path, organism: str, **_
) -> FlagEntry:
    REQUIRED_ANNOTATION_KEYS = {
        "SYMBOL",
        "GENENAME",
        "REFSEQ",
        "ENTREZID",
        "STRING_id",
        "GOSLIM_IDS",
        "EXTRA",
    }
    MASTER_ANNOTATION_KEY = {"_DEFAULT": "ENSEMBL", "Arabidopsis": "TAIR"}

    df_dge = pd.read_csv(dge_table)

    required_columns = REQUIRED_ANNOTATION_KEYS.union(
        {MASTER_ANNOTATION_KEY.get(organism, MASTER_ANNOTATION_KEY["_DEFAULT"])}
    )

    missing_columns = required_columns - set(df_dge.columns)
    # check logic
    if not missing_columns:
        code = FlagCode.GREEN
        message = f"Found all required annotation columns: {required_columns}"
    else:
        code = FlagCode.HALT
        message = (
            f"Missing the following required annotation columns: {missing_columns}"
        )
    return {"code": code, "message": message}


def check_dge_table_sample_columns_exist(
    dge_table: Path, samples: set[str], **_
) -> FlagEntry:
    # data specific preprocess
    df_dge = pd.read_csv(dge_table)

    missing_sample_columns = samples - set(df_dge.columns)

    # check logic
    if not missing_sample_columns:
        code = FlagCode.GREEN
        message = f"All samplewise columns present"
    else:
        code = FlagCode.HALT
        message = f"Missing these sample count columns: {missing_sample_columns}"
    return {"code": code, "message": message}


def check_dge_table_sample_columns_constraints(
    dge_table: Path, samples: set[str], **_
) -> FlagEntry:
    MINIMUM_COUNT = 0
    # data specific preprocess
    df_dge = pd.read_csv(dge_table)[samples]

    column_meets_constraints = df_dge.apply(
        lambda col: all(col >= MINIMUM_COUNT), axis="rows"
    )

    # check logic
    contraint_description = f"All counts are greater or equal to {MINIMUM_COUNT}"
    if all(column_meets_constraints):
        code = FlagCode.GREEN
        message = (
            f"All values in columns: {samples} met constraint: {contraint_description}"
        )
    else:
        code = FlagCode.HALT
        message = (
            f"These columns {list(column_meets_constraints.index[~column_meets_constraints])} "
            f"fail the contraint: {contraint_description}."
        )
    return {"code": code, "message": message}


def check_dge_table_group_columns_exist(
    dge_table: Path, runsheet: Path, **_
) -> FlagEntry:
    # data specific preprocess
    GROUP_PREFIXES = ["Group.Stdev_", "Group.Mean_"]
    expected_groups = utils_runsheet_to_expected_groups(runsheet)
    expected_columns = {
        "".join(comb)
        for comb in itertools.product(GROUP_PREFIXES, expected_groups.values())
    }
    df_dge_columns = set(pd.read_csv(dge_table).columns)
    missing_cols = expected_columns - df_dge_columns

    # check logic
    if not missing_cols:
        code = FlagCode.GREEN
        message = (
            f"All group summary stat columns present. {sorted(list(expected_columns))}"
        )
    else:
        code = FlagCode.HALT
        message = (
            f"Missing these group summary stat columns: {sorted(list(missing_cols))}"
        )
    return {"code": code, "message": message}


def check_dge_table_group_columns_constraints(
    dge_table: Path, runsheet: Path, **_
) -> FlagEntry:
    FLOAT_TOLERANCE = (
        0.001  # Percent allowed difference due to float precision differences
    )
    # data specific preprocess
    GROUP_PREFIXES = ["Group.Stdev_", "Group.Mean_"]
    expected_groups = utils_runsheet_to_expected_groups(runsheet)
    query_columns = {
        "".join(comb)
        for comb in itertools.product(GROUP_PREFIXES, expected_groups.values())
    }

    expected_group_lists = utils_runsheet_to_expected_groups(
        runsheet, map_to_lists=True
    )
    df_dge = pd.read_csv(dge_table)

    # issue trackers
    issues: dict[str, list[str]] = {
        f"mean computation deviates by more than {FLOAT_TOLERANCE} percent": [],
        f"standard deviation deviates by more than {FLOAT_TOLERANCE} percent": [],
    }

    group: str
    samples: list[str]
    for group, samples in expected_group_lists.items():
        abs_percent_differences = abs(
            (df_dge[f"Group.Mean_{group}"] - df_dge[samples].mean(axis="columns"))
            / df_dge[samples].mean(axis="columns")
            * 100
        )
        if any(abs_percent_differences > FLOAT_TOLERANCE):
            issues["mean computation incorrect"].append(group)

        abs_percent_differences = abs(
            (df_dge[f"Group.Stdev_{group}"] - df_dge[samples].std(axis="columns"))
            / df_dge[samples].mean(axis="columns")
            * 100
        )
        if any(abs_percent_differences > FLOAT_TOLERANCE):
            issues["standard deviation incorrect"].append(group)

    # check logic
    contraint_description = f"Group mean and standard deviations are correctly computed from samplewise normalized counts within a tolerance of {FLOAT_TOLERANCE} percent (to accomodate minor float related differences )"
    if not any([issue_type for issue_type in issues.values()]):
        code = FlagCode.GREEN
        message = f"All values in columns: {query_columns} met constraint: {contraint_description}"
    else:
        code = FlagCode.HALT
        message = (
            f"Issues found {issues} that"
            f"fail the contraint: {contraint_description}."
        )
    return {"code": code, "message": message}


def check_dge_table_comparison_statistical_columns_exist(
    dge_table: Path, runsheet: Path, **_
) -> FlagEntry:
    # data specific preprocess
    COMPARISON_PREFIXES = ["Log2fc_", "Stat_", "P.value_", "Adj.p.value_"]
    expected_groups = utils_runsheet_to_expected_groups(runsheet, map_to_lists=True)
    expected_comparisons = [
        "v".join(paired_groups)
        for paired_groups in itertools.permutations(expected_groups, 2)
    ]
    expected_columns = {
        "".join(comb)
        for comb in itertools.product(COMPARISON_PREFIXES, expected_comparisons)
    }
    df_dge_columns = set(pd.read_csv(dge_table).columns)
    missing_cols = expected_columns - df_dge_columns

    # check logic
    if not missing_cols:
        code = FlagCode.GREEN
        message = f"All comparision summary stat columns present. {sorted(list(expected_columns))}"
    else:
        code = FlagCode.HALT
        message = f"Missing these comparision summary stat columns: {sorted(list(missing_cols))}"
    return {"code": code, "message": message}


def utils_common_constraints_on_dataframe(
    df: pd.DataFrame, constraints: tuple[tuple[set, dict], ...]
) -> dict:

    issues: dict[str, list[str]] = {
        "Failed non null constraint": list(),
        "Failed non negative constraint": list(),
    }

    for (col_set, col_constraints) in constraints:
        # this will avoid overriding the original constraints dictionary
        # which is likely used in the check message
        col_constraints = col_constraints.copy()

        # limit to only columns of interest
        query_df = df[col_set]
        for (colname, colseries) in query_df.iteritems():
            # check non null constraint
            if col_constraints.pop("nonNull", False) and nonNull(colseries) == False:
                issues["Failed non null constraint"].append(colname)
            # check non negative constraint
            if (
                col_constraints.pop("nonNegative", False)
                and nonNegative(colseries) == False
            ):
                issues["Failed non negative constraint"].append(colname)
            # check allowed values constraint
            if allowedValues := col_constraints.pop("allowedValues", False):
                if onlyAllowedValues(colseries, allowedValues) == False:
                    issues["Failed non negative constraint"].append(colname)

            # raise exception if there are unhandled constraint keys
            if col_constraints:
                raise ValueError(f"Unhandled constraint types: {col_constraints}")

    return issues


def check_dge_table_group_statistical_columns_constraints(
    dge_table: Path, runsheet: Path, **_
) -> FlagEntry:
    expected_groups = utils_runsheet_to_expected_groups(runsheet, map_to_lists=True)
    expected_comparisons = [
        "v".join(paired_groups)
        for paired_groups in itertools.permutations(expected_groups, 2)
    ]

    resolved_constraints = (
        ({f"Log2fc_{comp}" for comp in expected_comparisons}, {"nonNull": True}),
        ({f"Stat_{comp}" for comp in expected_comparisons}, {"nonNull": True}),
        # can be removed from analysis before p-value and adj-p-value assessed
        # ref: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-are-some-p-values-set-to-na
        (
            {f"P.value_{comp}" for comp in expected_comparisons},
            {"nonNegative": True, "nonNull": False},
        ),
        (
            {f"Adj.p.value_{comp}" for comp in expected_comparisons},
            {"nonNegative": True, "nonNull": False},
        ),
    )

    df_dge = pd.read_csv(dge_table)

    # issue trackers
    # here: {prefix+constraint: [failed_columns]}
    issues: dict[str, list[str]] = dict()

    issues = utils_common_constraints_on_dataframe(df_dge, resolved_constraints)

    # check logic
    if not any([issue_type for issue_type in issues.values()]):
        code = FlagCode.GREEN
        message = f"All values in columns met constraint: {resolved_constraints}"
    else:
        code = FlagCode.HALT
        message = (
            f"Issues found {issues} that" f"fail the contraint: {resolved_constraints}."
        )
    return {"code": code, "message": message}


def check_dge_table_fixed_statistical_columns_exist(dge_table: Path, **_) -> FlagEntry:
    # data specific preprocess
    fixed_stats_columns = {
        "All.mean": {"nonNull": True, "nonNegative": True},
        "All.stdev": {"nonNull": True, "nonNegative": True},
        "LRT.p.value": {"nonNull": False, "nonNegative": True},
    }
    expected_columns = set(fixed_stats_columns)
    df_dge_columns = set(pd.read_csv(dge_table).columns)
    missing_cols = expected_columns - df_dge_columns

    # check logic
    if not missing_cols:
        code = FlagCode.GREEN
        message = f"All dataset summary stat columns present. {sorted(list(expected_columns))}"
    else:
        code = FlagCode.HALT
        message = (
            f"Missing these dataset summary stat columns: {sorted(list(missing_cols))}"
        )
    return {"code": code, "message": message}


def check_dge_table_fixed_statistical_columns_constraints(
    dge_table: Path, **_
) -> FlagEntry:
    # data specific preprocess
    fixed_stats_columns = (
        (["All.mean", "All.stdev"], {"nonNull": True, "nonNegative": True}),
        (["LRT.p.value"], {"nonNull": False, "nonNegative": True}),
    )

    df_dge = pd.read_csv(dge_table)

    # issue trackers
    # here: {prefix+constraint: [failed_columns]}
    issues: dict[str, list[str]] = dict()

    issues = utils_common_constraints_on_dataframe(df_dge, fixed_stats_columns)

    # check logic
    if not any([issue_type for issue_type in issues.values()]):
        code = FlagCode.GREEN
        message = f"All values in columns met constraint: {fixed_stats_columns}"
    else:
        code = FlagCode.HALT
        message = (
            f"Issues found {issues} that" f"fail the contraint: {fixed_stats_columns}."
        )
    return {"code": code, "message": message}


def check_dge_table_log2fc_within_reason(
    dge_table: Path, runsheet: Path, **_
) -> FlagEntry:
    LOG2FC_CROSS_METHOD_PERCENT_DIFFERENCE_THRESHOLD = 10  # Percent
    LOG2FC_CROSS_METHOD_TOLERANCE_PERCENT = 50  # Percent
    # data specific preprocess
    expected_groups = utils_runsheet_to_expected_groups(runsheet, map_to_lists=True)
    expected_comparisons = [
        "v".join(paired_groups)
        for paired_groups in itertools.permutations(expected_groups, 2)
    ]
    df_dge = pd.read_csv(dge_table)

    # Track error messages
    err_msg = ""
    for comparision in expected_comparisons:
        query_column = f"Log2fc_{comparision}"
        group1_mean_col = (
            "Group.Mean_" + comparision.split(")v(")[0] + ")"
        )  # Uses parens and adds them back to prevent slicing on 'v' within factor names
        group2_mean_col = "Group.Mean_" + "(" + comparision.split(")v(")[1]
        computed_log2fc = (df_dge[group1_mean_col] / df_dge[group2_mean_col]).apply(
            math.log, args=[2]
        )
        abs_percent_difference = abs(
            ((computed_log2fc - df_dge[query_column]) / df_dge[query_column]) * 100
        )
        percent_within_tolerance = (
            mean(
                abs_percent_difference
                < LOG2FC_CROSS_METHOD_PERCENT_DIFFERENCE_THRESHOLD
            )
            * 100
        )
        # flag if not enough within tolerance
        if percent_within_tolerance < LOG2FC_CROSS_METHOD_TOLERANCE_PERCENT:
            err_msg += (
                f"For comparison: '{comparision}' {percent_within_tolerance:.2f} % of genes have absolute percent differences "
                f"(between log2fc direct computation and DESeq2's approach) "
                f"less than {LOG2FC_CROSS_METHOD_PERCENT_DIFFERENCE_THRESHOLD} % which does not met the minimum percentage "
                f"({LOG2FC_CROSS_METHOD_TOLERANCE_PERCENT} %) of genes required.  "
                f"This may indicate misassigned or misaligned columns. "
            )

    if err_msg:
        code = FlagCode.YELLOW
        message = f"Pairwise statistical column(s) failing constraints: {err_msg}"
    else:
        code = FlagCode.GREEN
        message = (
            f"All log2fc within reason, specifically no more than {LOG2FC_CROSS_METHOD_TOLERANCE_PERCENT}% "
            f"of genes (actual %: {100 - percent_within_tolerance:.2f}) have a percent difference greater than "
            f"{LOG2FC_CROSS_METHOD_PERCENT_DIFFERENCE_THRESHOLD}%. "
        )

    return {"code": code, "message": message}


def check_viz_table_columns_exist(dge_table: Path, runsheet: Path, **_) -> FlagEntry:
    # data specific preprocess
    expected_groups = utils_runsheet_to_expected_groups(runsheet, map_to_lists=True)
    expected_comparisons = [
        "v".join(paired_groups)
        for paired_groups in itertools.permutations(expected_groups, 2)
    ]
    viz_pairwise_columns_prefixes = (
        (
            {f"Log2_Adj.p.value_{comp}" for comp in expected_comparisons},
            {"nonNull": False},
        ),
        (
            {f"Sig.1_{comp}" for comp in expected_comparisons},
            {"allowedValues": [False, True], "nonNull": False},
        ),
        (
            {f"Sig.05_{comp}" for comp in expected_comparisons},
            {"allowedValues": [False, True], "nonNull": False},
        ),
        (
            {f"Log2_P.value_{comp}" for comp in expected_comparisons},
            {"nonNegative": False, "nonNull": False},
        ),
        (
            {f"Updown_{comp}" for comp in expected_comparisons},
            {"allowedValues": [1, 0, -1], "nonNull": True},
        ),
    )

    expected_columns = set(
        itertools.chain(*[c1 for c1, _ in viz_pairwise_columns_prefixes])
    )
    df_dge_columns = set(pd.read_csv(dge_table).columns)
    missing_cols = expected_columns - df_dge_columns

    # check logic
    if not missing_cols:
        code = FlagCode.GREEN
        message = f"All viz specific comparison columns present. {sorted(list(expected_columns))}"
    else:
        code = FlagCode.HALT
        message = f"Missing these viz specific comparison columns: {sorted(list(missing_cols))}"
    return {"code": code, "message": message}


def check_viz_table_columns_constraints(
    dge_table: Path, runsheet: Path, **_
) -> FlagEntry:
    # data specific preprocess
    expected_groups = utils_runsheet_to_expected_groups(runsheet, map_to_lists=True)
    expected_comparisons = [
        "v".join(paired_groups)
        for paired_groups in itertools.permutations(expected_groups, 2)
    ]
    viz_pairwise_columns_constraints = (
        (
            {f"Log2_Adj.p.value_{comp}" for comp in expected_comparisons},
            {"nonNull": False},
        ),
        (
            {f"Sig.1_{comp}" for comp in expected_comparisons},
            {"allowedValues": [False, True], "nonNull": False},
        ),
        (
            {f"Sig.05_{comp}" for comp in expected_comparisons},
            {"allowedValues": [False, True], "nonNull": False},
        ),
        (
            {f"Log2_P.value_{comp}" for comp in expected_comparisons},
            {"nonNegative": False, "nonNull": False},
        ),
        (
            {f"Updown_{comp}" for comp in expected_comparisons},
            {"allowedValues": [1, 0, -1], "nonNull": True},
        ),
    )

    df_viz = pd.read_csv(dge_table)

    # issue trackers
    # here: {prefix+constraint: [failed_columns]}
    issues: dict[str, list[str]] = dict()

    issues = utils_common_constraints_on_dataframe(
        df_viz, viz_pairwise_columns_constraints
    )

    # check logic
    if not any([issue_type for issue_type in issues.values()]):
        code = FlagCode.GREEN
        message = (
            f"All values in columns met constraint: {viz_pairwise_columns_constraints}"
        )
    else:
        code = FlagCode.HALT
        message = (
            f"Issues found {issues} that"
            f"fail the contraint: {viz_pairwise_columns_constraints}."
        )
    return {"code": code, "message": message}


def check_viz_pca_table_index_and_columns_exist(
    pca_table: Path, samples: set[str]
) -> FlagEntry:
    EXPECTED_VIS_PCA_COLUMNS = {"PC1", "PC2", "PC3"}
    err_msg = ""
    # data specific preprocess
    df = pd.read_csv(pca_table, index_col=0)

    # check all samples included
    if missing_samples := samples - set(df.index):
        err_msg += f"Missing samples in index: {missing_samples}"

    # check all expected columns exist
    if missing_cols := EXPECTED_VIS_PCA_COLUMNS - set(df.columns):
        err_msg += f"Missing expected columns: {missing_cols}"

    if not err_msg:
        code = FlagCode.GREEN
        message = f"PCA Table has all the samples in the index and these columns exist: {EXPECTED_VIS_PCA_COLUMNS}"
    else:
        code = FlagCode.HALT
        message = err_msg

    return {"code": code, "message": message}

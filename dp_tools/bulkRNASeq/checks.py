from collections import defaultdict
import copy
import enum
import gzip
import logging
import math
from pathlib import Path
from statistics import mean, median, stdev
import subprocess
from typing import Callable, DefaultDict, Dict, List, Set, Tuple, Union
import numpy as np

import pandas as pd
from dp_tools.components.components import GenomeAlignments, RawReadsComponent
from dp_tools.core.entity_model import (
    BaseComponent,
    DataDir,
    DataFile,
    ModuleLevelMQC,
    TemplateDataset,
    TemplateSample,
)

log = logging.getLogger(__name__)

from dp_tools.core.check_model import Check, Flag, FlagCode, FlagEntry, Flaggable

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


class MIDDLE(enum.Enum):
    mean: Tuple[Callable] = (mean,)
    median: Tuple[Callable] = (median,)

    def __call__(self, *args, **kwargs):
        return self.value[0](*args, **kwargs)


def identify_outliers(
    valueDict: Dict[str, float], standard_deviation_threshold: float, middle: Callable
):
    # determine middle value
    middle_value: float = middle(valueDict.values())
    std_deviation: float = stdev(valueDict.values())

    # init tracker
    # holds the key name and the standard deviations from the middle
    outliers: Dict[str, float] = dict()

    # exit early if std_deviation is zero (i.e. no outliers)
    if std_deviation == 0:
        return outliers

    # check if a value is an outlier
    for key, value in valueDict.items():
        # calculate standard deviations
        num_std_deviations_vector = (value - middle_value) / std_deviation
        # if an outlier, add it to a dict of outliers (include a +/- standard deviations)
        if abs(num_std_deviations_vector) > standard_deviation_threshold:
            outliers[key] = num_std_deviations_vector

    return outliers


# TODO: typedict for thresholds
def identify_values_past_thresholds(thresholds: dict, value: float) -> List[FlagCode]:
    """Return empty list if no codes are raised"""
    VALID_THRESHOLD_TYPES = {"lower", "upper"}
    new_codes = list()
    for threshold in thresholds:
        assert (
            threshold.get("type") in VALID_THRESHOLD_TYPES
        ), f"Invalid threshold type configured: valid options {VALID_THRESHOLD_TYPES} got {threshold.get('type')}"
        if threshold.get("type") == "lower":
            if value < threshold["value"]:
                new_codes.append(threshold["code"])
        elif threshold.get("type") == "upper":
            if value > threshold["value"]:
                new_codes.append(threshold["code"])
    return new_codes


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


class DATASET_RAWREADS_0001(Check):
    config = {
        "metrics": [
            "percent_gc",
            "avg_sequence_length",
            "total_sequences",
            "percent_duplicates",
            # "percent_fails", number of failed FastQC submodules, not a very useful metric for BulkRNASeq
        ],
        "middle": MIDDLE.median,
        "yellow_standard_deviation_threshold": 2,
        "red_standard_deviation_threshold": 4,
        "target_components_by_paired_end": {
            True: ["rawForwardReads", "rawReverseReads"],
            False: ["rawReads"],
        },
    }
    description = (
        "Check that the reads stats (source from FastQC) have no outliers among samples "
        "for the following metrics: {metrics}. "
        "Yellow Flagged Outliers are defined as a being {yellow_standard_deviation_threshold} - {red_standard_deviation_threshold} standard "
        "deviations away from the {middle.name}. "
        "Red Flagged Outliers are defined as a being {red_standard_deviation_threshold}+ standard "
        "deviations away from the {middle.name}. "
    )
    flag_desc = {
        FlagCode.GREEN: "No reads metric outliers detected for {metrics}",
        FlagCode.YELLOW1: "Outliers detected as follows (values are rounded number of standard deviations from middle): {formatted_outliers}",
        FlagCode.RED1: "Outliers detected as follows (values are rounded number of standard deviations from middle): {formatted_outliers}",
    }

    def validate_func(self: Check, dataset: TemplateDataset) -> Flag:
        code = FlagCode.GREEN

        # pull variables from config
        metrics = self.config["metrics"]
        middle = self.config["middle"]
        yellow_threshold = self.config["yellow_standard_deviation_threshold"]
        red_threshold = self.config["red_standard_deviation_threshold"]

        # init trackers for issues
        outliers: DefaultDict[str, Dict[str, float]] = defaultdict(dict)

        # determine reads components in samples
        readsComponents = self.config["target_components_by_paired_end"][
            dataset.metadata.paired_end
        ]

        def format_identifier(sample_name: str, component_str: str) -> str:
            """Add forward and reverse suffix if paired end, add nothing otherwise"""
            return (
                f"{sample_name}:{component_str}"
                if dataset.metadata.paired_end
                else sample_name
            )

        # iterate through metrics (here all pulled from FastQC general stats)
        for readComponent in readsComponents:
            for metric in metrics:
                sampleToMetric: Dict[str, float] = {
                    format_identifier(s.name, readComponent): getattr(
                        s, readComponent
                    ).mqcData["FastQC"]["General_Stats"][metric]
                    for s in dataset.samples.values()
                }

                # ensure any NaN convert to zero as implied by MultiQC
                sampleToMetric = convert_nan_to_zero(sampleToMetric)

                # yellow level outliers
                if outliersForThisMetric := identify_outliers(
                    sampleToMetric,
                    standard_deviation_threshold=yellow_threshold,
                    middle=middle,
                ):
                    if code < FlagCode.YELLOW1:
                        code = FlagCode.YELLOW1
                    outliers[metric] = outliers[metric] | outliersForThisMetric

                # red level outliers
                if outliersForThisMetric := identify_outliers(
                    sampleToMetric,
                    standard_deviation_threshold=red_threshold,
                    middle=middle,
                ):
                    if code < FlagCode.RED1:
                        code = FlagCode.RED1
                    outliers[metric] = outliers[metric] | outliersForThisMetric

        return Flag(
            codes=code,
            check=self,
            message_args={
                "outliers": outliers,
                "metrics": metrics,
                "formatted_outliers": pformat(outliers, formatfloat),
            },
        )


class DATASET_TRIMREADS_0001(DATASET_RAWREADS_0001):
    # overwrite specific config only
    config = DATASET_RAWREADS_0001.config | {
        "target_components_by_paired_end": {
            True: ["trimForwardReads", "trimReverseReads"],
            False: ["trimReads"],
        }
    }


class DATASET_GENOMEALIGNMENTS_0001(Check):
    config = {
        "metrics": [
            # "total_reads", # check in FastQC, but is used to normalize
            # "avg_input_read_length",
            # "uniquely_mapped", # redundant with better metric of percent
            "uniquely_mapped_percent",
            "avg_mapped_read_length",
            # "num_splices",
            # "num_annotated_splices",
            # "num_GTAG_splices",
            # "num_GCAG_splices",
            # "num_ATAC_splices",
            # "num_noncanonical_splices",
            "mismatch_rate",
            "deletion_rate",
            "deletion_length",
            "insertion_rate",
            "insertion_length",
            # "multimapped", # redundant with better metric of percent
            "multimapped_percent",
            # "multimapped_toomany",  # redundant with better metric of percent
            "multimapped_toomany_percent",
            "unmapped_mismatches_percent",
            "unmapped_tooshort_percent",
            "unmapped_other_percent",
            # "unmapped_mismatches", # redundant with better metric of percent
            # "unmapped_tooshort", # redundant with better metric of percent
            # "unmapped_other", # redundant with better metric of percent
        ],
        "middle": MIDDLE.median,
        "yellow_standard_deviation_threshold": 2,
        "red_standard_deviation_threshold": 4,
    }
    description = (
        "Check that the genome alignment stats (source from STAR logs) have no outliers among samples "
        "for the following metrics: {metrics}. "
        "Yellow Flagged Outliers are defined as a being {yellow_standard_deviation_threshold} - {red_standard_deviation_threshold} standard "
        "deviations away from the {middle.name}. "
        "Red Flagged Outliers are defined as a being {red_standard_deviation_threshold}+ standard "
        "deviations away from the {middle.name}. "
    )
    flag_desc = {
        FlagCode.GREEN: "No genome alignment metric outliers detected for {metrics}",
        FlagCode.YELLOW1: "Outliers detected as follows (values are rounded number of standard deviations from middle): {formatted_outliers}",
        FlagCode.RED1: "Outliers detected as follows (values are rounded number of standard deviations from middle): {formatted_outliers}",
    }

    def validate_func(self: Check, dataset: TemplateDataset) -> Flag:
        code = FlagCode.GREEN

        # pull variables from config
        metrics = self.config["metrics"]
        middle = self.config["middle"]
        yellow_threshold = self.config["yellow_standard_deviation_threshold"]
        red_threshold = self.config["red_standard_deviation_threshold"]

        # init trackers for issues
        outliers: DefaultDict[str, Dict[str, float]] = defaultdict(dict)

        # determine reads components in samples
        targetComponents = ["genomeAlignments"]

        # iterate through metrics (here all pulled from FastQC general stats)
        for targetComponent in targetComponents:
            for metric in metrics:
                sampleToMetric: Dict[str, float] = {
                    s.name: getattr(s, targetComponent).mqcData["STAR"][
                        "General_Stats"
                    ][metric]
                    for s in dataset.samples.values()
                }
                # ensure any NaN convert to zero as implied by MultiQC
                sampleToMetric = convert_nan_to_zero(sampleToMetric)

                # yellow level outliers
                if outliersForThisMetric := identify_outliers(
                    sampleToMetric,
                    standard_deviation_threshold=yellow_threshold,
                    middle=middle,
                ):
                    if code < FlagCode.YELLOW1:
                        code = FlagCode.YELLOW1
                    outliers[metric] = outliers[metric] | outliersForThisMetric

                # red level outliers
                if outliersForThisMetric := identify_outliers(
                    sampleToMetric,
                    standard_deviation_threshold=red_threshold,
                    middle=middle,
                ):
                    if code < FlagCode.RED1:
                        code = FlagCode.RED1
                    outliers[metric] = outliers[metric] | outliersForThisMetric

        return Flag(
            codes=code,
            check=self,
            message_args={
                "outliers": outliers,
                "metrics": metrics,
                "formatted_outliers": pformat(outliers, formatfloat),
            },
        )


class DATASET_RSEQCANALYSIS_0001(Check):
    config = {
        "plots_all": ["Read Distribution", "Infer experiment", "Gene Body Coverage"],
        "plot_paired_end": ["Inner Distance"],
        "middle": MIDDLE.median,
        "yellow_standard_deviation_threshold": 2,
        "red_standard_deviation_threshold": 4,
        "stranded_assessment_range": {"min": 75, "max": 100},  # percents
        "halt_ambiguous_dominant_strandedness_range": {
            "min": 60,
            "max": 75,
        },  # percents
        "unstranded_assessment_range": {"min": 40, "max": 60},  # percents
        "valid_dominant_strandedness_assessments": [
            "Sense (% Tags)",
            "Antisense (% Tags)",
        ],  # this leaves out undetermined, which should raise alarms if it is the dominant assessment
    }
    description = (
        "Check that the rseqc analysis stats (sourced from the rseqc logs) have no outlier values among samples "
        "for the following plots: {plots_all} (Paired end only: {plot_paired_end}). "
        "Yellow Flagged Outliers are defined as a being {yellow_standard_deviation_threshold} - {red_standard_deviation_threshold} standard "
        "deviations away from the {middle.name}. "
        "Red Flagged Outliers are defined as a being {red_standard_deviation_threshold}+ standard "
        "deviations away from the {middle.name}. "
        "Additionally the following is assessed for infer experiment strandedess metrics: "
        "A Halt Flag is raised in the case that the dominant strandessness is between "
        "{halt_ambiguous_dominant_strandedness_range} "
        "Note: the 'dominant strandedness' is the max(datasetwide_median(antisense), datasetwide_median(sense)) "
        "Valid assessments include {valid_dominant_strandedness_assessments}, other assessments (e.g. 'undetermined') will raise a Halting flag "
    )
    flag_desc = {
        FlagCode.GREEN: "No rseqc analysis metric outliers detected for {metrics}",
        FlagCode.YELLOW1: "Outliers detected as follows (values are rounded number of standard deviations from middle): {formatted_outliers}",
        FlagCode.RED1: "Outliers detected as follows (values are rounded number of standard deviations from middle): {formatted_outliers}",
        FlagCode.RED2: "At least one sample is outside the dominant strandedness assignment range: {samples_outside_range}",
        FlagCode.HALT1: "The dominant strandedness is {dominant_strandedness}, this is lower than the halting flag threshold.",
        FlagCode.HALT2: "The dominant strandedness is {dominant_strandedness} which is not a invalid assessment.",
    }

    def validate_func(self: Check, dataset: TemplateDataset) -> Flag:
        codes = {FlagCode.GREEN}

        # pull variables from config
        targetPlotsAll = self.config["plots_all"]
        targetPlotsPairedEnd = self.config["plot_paired_end"]
        middle = self.config["middle"]
        yellow_threshold = self.config["yellow_standard_deviation_threshold"]
        red_threshold = self.config["red_standard_deviation_threshold"]

        # init trackers for issues
        outliers: DefaultDict[str, Dict[str, float]] = defaultdict(dict)

        # extend with paired end specific plot if appropriate
        targetPlots = targetPlotsAll
        if dataset.metadata.paired_end:
            targetPlots.extend(targetPlotsPairedEnd)

        # iterate through metrics (here all pulled from FastQC general stats)
        for plot_name in targetPlots:
            # extract dataframe of all samples
            df = dataset.getMQCDataFrame(
                sample_component="rSeQCAnalysis", mqc_module="RSeQC", mqc_plot=plot_name
            )

            # convert to samplewise dicts
            metricToSampleToMetricValue: Dict[str, Dict[str, float]] = df.to_dict()

            for metricName, sampleToMetricValue in metricToSampleToMetricValue.items():
                # ensure any NaN convert to zero as implied by MultiQC
                sampleToMetricValue = convert_nan_to_zero(sampleToMetricValue)

                # yellow level outliers
                if outliersForThisMetric := identify_outliers(
                    sampleToMetricValue,
                    standard_deviation_threshold=yellow_threshold,
                    middle=middle,
                ):
                    if max(codes) < FlagCode.YELLOW1:
                        codes.add(FlagCode.YELLOW1)
                    outliers[metricName] = outliers[metricName] | outliersForThisMetric

                # red level outliers
                if outliersForThisMetric := identify_outliers(
                    sampleToMetricValue,
                    standard_deviation_threshold=red_threshold,
                    middle=middle,
                ):
                    if max(codes) < FlagCode.RED1:
                        codes.add(FlagCode.RED1)
                        # remove lower FlagCode YELLOW1
                        codes.remove(FlagCode.YELLOW1)
                    outliers[metricName] = outliers[metricName] | outliersForThisMetric

        def get_median_strandedness(dataset: TemplateDataset) -> tuple[str, float]:
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
        if (
            strand_assessment
            not in self.config["valid_dominant_strandedness_assessments"]
        ):
            codes.add(FlagCode.HALT2)

        # flag based on thresholds
        assessment_value: float = median_strandedness[strand_assessment]

        is_stranded: bool = (
            self.config["stranded_assessment_range"]["max"]
            > assessment_value
            > self.config["stranded_assessment_range"]["min"]
        )
        is_unstranded: bool = (
            self.config["unstranded_assessment_range"]["max"]
            > assessment_value
            > self.config["unstranded_assessment_range"]["min"]
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
                self.config["stranded_assessment_range"]["min"],
                self.config["stranded_assessment_range"]["max"],
            )
        elif is_unstranded:
            samples_outside_range = determine_samples_outside_range(
                dataset,
                self.config["unstranded_assessment_range"]["min"],
                self.config["unstranded_assessment_range"]["max"],
            )
        else:  # this means that the standing is ambiguous
            samples_outside_range = list()
            codes.add(FlagCode.HALT1)

        if len(samples_outside_range) != 0:
            codes.add(FlagCode.RED2)

        return Flag(
            codes=codes,
            check=self,
            message_args={
                "outliers": outliers,
                "formatted_outliers": pformat(outliers, formatfloat),
                "dominant_strandedness": (strand_assessment, assessment_value),
                "samples_outside_range": samples_outside_range,
            },
        )


class DATASET_GENECOUNTS_0001(Check):
    config = {
        "metrics": [
            "Unalignable",
            "Alignable",
            "Filtered",
            "Total",
            "alignable_percent",
            "Unique",
            "Multi",
            "Uncertain",
        ],
        "middle": MIDDLE.median,
        "yellow_standard_deviation_threshold": 2,
        "red_standard_deviation_threshold": 4,
    }
    description = (
        "Check that the gene counts alignments (source from the RSEM logs) have no outlier values among samples "
        "for the following metrics: {metrics} "
        "Yellow Flagged Outliers are defined as a being {yellow_standard_deviation_threshold} - {red_standard_deviation_threshold} standard "
        "deviations away from the {middle.name}. "
        "Red Flagged Outliers are defined as a being {red_standard_deviation_threshold}+ standard "
        "deviations away from the {middle.name}. "
    )
    flag_desc = {
        FlagCode.GREEN: "No gene count mapping metric outliers detected for {metrics}",
        FlagCode.YELLOW1: "Outliers detected as follows (values are rounded number of standard deviations from middle): {formatted_outliers}",
        FlagCode.RED1: "Outliers detected as follows (values are rounded number of standard deviations from middle): {formatted_outliers}",
    }

    def validate_func(self: Check, dataset: TemplateDataset) -> Flag:
        codes = {FlagCode.GREEN}

        # pull variables from config
        metrics = self.config["metrics"]
        middle = self.config["middle"]
        yellow_threshold = self.config["yellow_standard_deviation_threshold"]
        red_threshold = self.config["red_standard_deviation_threshold"]

        # init trackers for issues
        outliers: DefaultDict[str, Dict[str, float]] = defaultdict(dict)

        # extract dataframe of general stats
        df = dataset.getMQCDataFrame(
            sample_component="geneCounts", mqc_module="Rsem", mqc_plot="general_stats"
        )

        # iterate through metrics (here all pulled from FastQC general stats)
        for metric_name in metrics:
            sampleToMetricValue = df[[metric_name]].to_dict()[metric_name]
            # ensure any NaN convert to zero as implied by MultiQC
            sampleToMetricValue = convert_nan_to_zero(sampleToMetricValue)

            # yellow level outliers
            if outliersForThisMetric := identify_outliers(
                sampleToMetricValue,
                standard_deviation_threshold=yellow_threshold,
                middle=middle,
            ):
                if max(codes) < FlagCode.YELLOW1:
                    codes.add(FlagCode.YELLOW1)
                    outliers[metric_name] = (
                        outliers[metric_name] | outliersForThisMetric
                    )

            # red level outliers
            if outliersForThisMetric := identify_outliers(
                sampleToMetricValue,
                standard_deviation_threshold=red_threshold,
                middle=middle,
            ):
                if max(codes) < FlagCode.RED1:
                    codes.add(FlagCode.RED1)
                    outliers[metric_name] = (
                        outliers[metric_name] | outliersForThisMetric
                    )

        return Flag(
            codes=codes,
            check=self,
            message_args={
                "outliers": outliers,
                "formatted_outliers": pformat(outliers, formatfloat),
                "metrics": metrics,
            },
        )


# TODO: Flag message gets really messy, convert into a json like string for easier reading/parsing
# TODO: Check for extra unexpected columns, these should give clues to names differences
class DATASET_DIFFERENTIALGENEEXPRESSION_CONTRASTS_0001(Check):
    USE_SUBCHECK_MONADS = True

    config = {
        "expected_tables": [
            "differential_expression.csv",
            "visualization_output_table.csv",
            "visualization_PCA_table.csv",
        ],
        # Expected column name, but dependent on dataset organism
        "dge_table_master_annotation_keys": {
            "Arabidopsis thaliana": "TAIR",
            "_DEFAULT": "ENSEMBL",
        },
        "dge_table_expected_annotation_columns": [
            "SYMBOL",
            "GENENAME",
            "REFSEQ",
            "ENTREZID",
            "STRING_id",
            "GOSLIM_IDS",
        ],
        # includes column specific constraints
        # these prefix as follows {prefix}{pairWiseFactorGroupComparison}
        "pairwise_columns_prefixes": {
            "Log2fc_": {"nonNull": True},
            "Stat_": {"nonNull": True},
            # can be removed from analysis before p-value and adj-p-value assessed
            # ref: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-are-some-p-values-set-to-na
            "P.value_": {"nonNegative": True, "nonNull": False},
            "Adj.p.value_": {"nonNegative": True, "nonNull": False},
        },
        "viz_pairwise_columns_prefixes": {
            "Log2_Adj.p.value_": {"nonNull": False},
            "Sig.1_": {"allowedValues": [False, True], "nonNull": False},
            "Sig.05_": {"allowedValues": [False, True], "nonNull": False},
            "Log2_P.value_": {"nonNegative": False, "nonNull": False},
            "Updown_": {"allowedValues": [1, 0, -1], "nonNull": True},
        },
        # these prefix as follows {prefix}{FactorGroup}
        "group_factorwise_columns_prefixes": {
            "Group.Mean_": {"nonNull": True, "nonNegative": True},
            "Group.Stdev_": {"nonNull": True, "nonNegative": True},
        },
        "fixed_stats_columns": {
            "All.mean": {"nonNull": True, "nonNegative": True},
            "All.stdev": {"nonNull": True, "nonNegative": True},
            "LRT.p.value": {"nonNull": False, "nonNegative": True},
        },
        "sample_counts_constraints": {"nonNegative": True},
        "expected_vis_pca_columns": [
            "PC1",
            "PC2",
        ],  # more may be included but these are REQUIRED
        "float_tolerance": 0.0001,  # PERCENT
        # TODO: DISCUSS, these baseline values, should indicate a very heavy left-hand skewed histogram of differences - JDO
        "log2fc_cross_method_percent_difference_threshold": 10,  # PERCENT
        "log2fc_cross_method_tolerance_percent": 50,  # PERCENT
        # PERCENT difference minimum between group means to included
        "log2fc_cross_method_sign_check_group_mean_difference_threshold": 50,
        # PERCENT genes allowed sign inversions between methods for groups that meet
        #   log2fc_cross_method_sign_check_group_mean_difference_threshold minimum
        "log2fc_cross_method_sign_check_tolerance_percent": 0,
        # "middle": MIDDLE.median,
        # "yellow_standard_deviation_threshold": 2,
        # "red_standard_deviation_threshold": 4,
    }
    description = (
        "Check that the differential expression outputs exist (source from the deseq2 script) and  "
        "the following tables: {expected_tables}.  "
        "For studies with ERCC spike-in, performs the same check on analogous tables. "
        "Additional performs the file specific validations: "
        "- contrasts.csv: Includes all the existing comparison groups (based on factor values in the metadata) and is formatted correctly"
        "- differential_expression.csv:  Includes expected annotation columns {dge_table_expected_annotation_columns}, includes a master annotation key "
        "column dependent on the dataset organism as follows: {dge_table_master_annotation_keys} ,"
        "includes sample count columns for all samples, all sample count values are non-negative, "
        "all pairwise comparision columns exist with the following prefixes and adhere to the following constraints: {pairwise_columns_prefixes} "
        "all groupFactorWise statistics columns exists with the following prefixes and adhere to the following constraints: {group_factorwise_columns_prefixes} "
        "all fixed statistics columns exist and adhere to the following constraints: {fixed_stats_columns} "
        " - visualization_PCA_table.csv: All samples in index and at the following columns exist {expected_vis_pca_columns} "
        " - visualization_output_table.csv: Performs same checks as differential_expression.csv as well as, "
        "ensuring the additional pairwise comparision columns exist with the following prefixes and "
        "adhere to the following constraints: {expected_vis_pca_columns}. "
        "Confirms that gene counts between differential expression table and normalized counts tables are the same. "
        "Confirms that computations match expectations with respect to following operations: (Float tolerance: +/-{float_tolerance} %)"
        "- Group means are correctly computed from normalized counts "
        "- log2FC values (computed with DESeq2's MLE approach) are comparable to direct computation with log2( mean(group1) / mean(group2) ), specifically "
        "checking if at least {log2fc_cross_method_tolerance_percent} % of all genes have absolute percent differences between methods "
        "less than {log2fc_cross_method_percent_difference_threshold} % "
    )
    flag_desc = {}

    def validate_func(self: Check, dataset: TemplateDataset) -> Flag:
        all_flags: list[FlagEntry] = list()

        contrasts_flag_data = (
            Flaggable(
                dataset=dataset,
                componentTarget="differentialGeneExpression",
            )
            .bind(self.check_contrasts_headers)
            .bind(self.check_contrasts_rows)
        ).export_flag_data_to_Flag(check=self)

        return contrasts_flag_data

        all_flags.extend(
            (
                Flaggable(
                    dataset=dataset,
                    componentTarget="differentialGeneExpression",
                )
                .bind(self.check_contrasts_headers)
                .bind(self.check_contrasts_rows)
            ).flag_data
        )
        # wrap inputs with flag_data
        all_flags.extend(
            (
                Flaggable(
                    dataset=dataset,
                    componentTarget="differentialGeneExpression",
                    componentDataAsset="annotatedTableCSV",
                )
                .bind(self.check_dge_table_annotation_columns_exist)
                .bind(self.check_dge_table_sample_columns_exist)
                .bind(self.check_dge_table_sample_column_constraints)
                .bind(self.check_dge_table_statistical_columns_exist)
                .bind(self.check_dge_table_statistical_column_constraints)
                .bind(self.check_dge_table_group_factor_columns_exist)
                .bind(self.check_dge_table_group_factor_column_constraints)
                .bind(self.check_dge_table_fixed_statistical_columns_exist)
                .bind(self.check_dge_table_fixed_statistical_column_constraints)
                .bind(self.check_dge_table_log2fc_within_reason)
                .bind(self.check_dge_table_group_means_computation)
            ).flag_data
        )
        all_flags.extend(
            (
                Flaggable(
                    dataset=dataset,
                    componentTarget="differentialGeneExpression",
                    componentDataAsset="visualizationTableCSV",
                )
                .bind(self.check_dge_table_annotation_columns_exist)
                .bind(self.check_dge_table_sample_columns_exist)
                .bind(self.check_dge_table_sample_column_constraints)
                .bind(self.check_dge_table_statistical_columns_exist)
                .bind(self.check_dge_table_statistical_column_constraints)
                .bind(self.check_dge_table_group_factor_columns_exist)
                .bind(self.check_dge_table_group_factor_column_constraints)
                .bind(self.check_dge_table_fixed_statistical_columns_exist)
                .bind(self.check_dge_table_fixed_statistical_column_constraints)
                .bind(self.check_dge_table_log2fc_within_reason)
                .bind(self.check_dge_table_group_means_computation)
                .bind(self.check_viz_pca_table_index_and_columns_exist)
            ).flag_data
        )
        if dataset.metadata.has_ercc:
            all_flags.extend(
                (
                    Flaggable(
                        dataset=dataset,
                        componentTarget="differentialGeneExpressionERCC",
                    )
                    .bind(self.check_contrasts_headers)
                    .bind(self.check_contrasts_rows)
                ).flag_data
            )
            all_flags.extend(
                (
                    Flaggable(
                        dataset=dataset,
                        componentTarget="differentialGeneExpressionERCC",
                        componentDataAsset="annotatedTableCSV",
                    )
                    .bind(self.check_dge_table_annotation_columns_exist)
                    .bind(self.check_dge_table_sample_columns_exist)
                    .bind(self.check_dge_table_sample_column_constraints)
                    .bind(self.check_dge_table_statistical_columns_exist)
                    .bind(self.check_dge_table_statistical_column_constraints)
                    .bind(self.check_dge_table_group_factor_columns_exist)
                    .bind(self.check_dge_table_group_factor_column_constraints)
                    .bind(self.check_dge_table_fixed_statistical_columns_exist)
                    .bind(self.check_dge_table_fixed_statistical_column_constraints)
                    .bind(self.check_dge_table_log2fc_within_reason)
                    .bind(self.check_dge_table_group_means_computation)
                ).flag_data
            )
            all_flags.extend(
                (
                    Flaggable(
                        dataset=dataset,
                        componentTarget="differentialGeneExpressionERCC",
                        componentDataAsset="visualizationTableCSV",
                    )
                    .bind(self.check_dge_table_annotation_columns_exist)
                    .bind(self.check_dge_table_sample_columns_exist)
                    .bind(self.check_dge_table_sample_column_constraints)
                    .bind(self.check_dge_table_statistical_columns_exist)
                    .bind(self.check_dge_table_statistical_column_constraints)
                    .bind(self.check_dge_table_group_factor_columns_exist)
                    .bind(self.check_dge_table_group_factor_column_constraints)
                    .bind(self.check_dge_table_fixed_statistical_columns_exist)
                    .bind(self.check_dge_table_fixed_statistical_column_constraints)
                    .bind(self.check_dge_table_log2fc_within_reason)
                    .bind(self.check_dge_table_group_means_computation)
                    .bind(self.check_viz_pca_table_index_and_columns_exist)
                ).flag_data
            )
        # unwrap and return flag_data
        print(1)

    @staticmethod
    def check_contrasts_headers(
        dataset: TemplateDataset, componentTarget: str
    ) -> FlagEntry:
        # extract target Component
        target_component = getattr(dataset, componentTarget)
        # extract dicts for deseq2 contrasts and the metadata formatted one here
        # make sure to read in explicit index column for deseq2
        dict_deseq2: Dict = pd.read_csv(
            target_component.contrastsCSV.path, index_col=0
        ).to_dict(orient="list")
        dict_data_model: Dict = dataset.metadata.contrasts.to_dict(orient="list")

        # check that all headers are present
        deseq2_headers = set(dict_deseq2.keys())
        data_model_headers = set(dict_data_model.keys())
        if deseq2_headers != data_model_headers:
            code = FlagCode.HALT1
            message = (
                f"Contrast headers do not match expecations "
                f"(based on Factor Values in {dataset.metadata.runsheet.path.name}). "
                f"Mismatched headers: (runsheet-unique: {data_model_headers - deseq2_headers}) "
                f"(contrasts-file-unique: {deseq2_headers - data_model_headers})"
            )
        else:
            code = FlagCode.GREEN
            message = f"Contrast headers correctly formed (based on Factor Values in {dataset.metadata.runsheet.path.name})"
        return {"code": code, "message": message}

    @staticmethod
    def check_contrasts_rows(
        dataset: TemplateDataset, componentTarget: str
    ) -> FlagEntry:
        # extract target Component
        target_component = getattr(dataset, componentTarget)
        # extract dicts for deseq2 contrasts and the metadata formatted one here
        # make sure to read in explicit index column for deseq2
        dict_deseq2: Dict = pd.read_csv(
            target_component.contrastsCSV.path, index_col=0
        ).to_dict(orient="list")
        dict_data_model: Dict = dataset.metadata.contrasts.to_dict(orient="list")

        # check contents of each column matches expecatation (group1 and group2 formatted as expected)
        # this also rechecks headers (keys) but that is caught in the prior validation
        if dict_deseq2 != dict_data_model:
            code = FlagCode.HALT1
            message = (
                f"Rows don't match expectations. Deseq2: {dict_deseq2}. "
                f"DataModel (from metadata source): {dict_data_model}"
            )
        else:
            code = FlagCode.GREEN
            message = f"Contrast rows correctly formed (based on Factor Values in {dataset.metadata.runsheet.path.name})"
        return {"code": code, "message": message}

    def check_dge_table_annotation_columns_exist(
        self,
        dataset: TemplateDataset,
        componentTarget: str,
        componentDataAsset: str = "annotatedTableCSV",
    ):
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)

        # read in dataframe
        df_dge = pd.read_csv(target_data_asset.path)

        # check all constant columns exist
        missing_constant_columns: set
        master_key = self.config["dge_table_master_annotation_keys"].get(
            dataset.metadata.organism,
            self.config["dge_table_master_annotation_keys"]["_DEFAULT"],
        )
        log.debug(
            f"Resolved master annotation key for {dataset.metadata.organism} is {master_key}"
        )
        expected_columns: list = self.config["dge_table_expected_annotation_columns"] + [master_key]  # type: ignore
        if missing_constant_columns := set(expected_columns) - set(df_dge.columns):
            code = FlagCode.HALT2
            err_msg += f"Annotation Columns missing: {missing_constant_columns}"
        else:
            code = FlagCode.GREEN
            message = f"All the following expected annotation columns exist: {expected_columns}"

        return {"code": code, "message": message}

    @staticmethod
    def check_dge_table_sample_columns_exist(
        dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)
        df_dge = pd.read_csv(target_data_asset.path)
        # check all sample counts columns exist
        expected_samples = set(dataset.samples.keys())
        if missing_samples := expected_samples - set(df_dge.columns):
            code = FlagCode.HALT2
            message = f"Sample Count Columns missing: {missing_samples}"
        else:
            code = FlagCode.GREEN
            message = f"Sample Count Columns exist as follows: {list(df_dge.columns)}"

        # check that they met constraints
        # all sample column counts are not negative
        if not (df_dge[list(expected_samples)] >= 0).all(axis=None):
            err_msg += (
                f"Sample Count Columns include negative values: {missing_samples}"
            )

        return {"code": code, "message": message}

    @staticmethod
    def check_dge_table_sample_column_constraints(
        dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)
        df_dge = pd.read_csv(target_data_asset.path)

        expected_samples = set(dataset.samples.keys())

        # check that they met constraints
        # all sample column counts are not negative
        if not (df_dge[list(expected_samples)] >= 0).all(axis=None):
            code = FlagCode.HALT2
            message = f"Sample Count Columns include negative values"
        else:
            code = FlagCode.GREEN
            message = f"Sample Count Columns meet constraints: 'Non-negative'"

        return {"code": code, "message": message}

    def check_dge_table_statistical_columns_exist(
        self, dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)
        df_dge = pd.read_csv(target_data_asset.path)

        aggregate_missing_cols: list = list()
        aggregate_target_cols: list = list()

        # check all expected statistic columns present
        # pairwise comparison level
        pairwise_comparisons = dataset.metadata.contrasts.columns
        for statistical_prefix in self.config[
            "pairwise_columns_prefixes"
        ].keys():  # type: ignore
            target_cols: list = [
                f"{statistical_prefix}{comparison}"
                for comparison in pairwise_comparisons
            ]
            aggregate_target_cols += target_cols
            # check existense first and bail if any don't exist
            if missing_cols := set(target_cols) - set(df_dge.columns):
                aggregate_missing_cols += missing_cols

        if aggregate_missing_cols:
            code = FlagCode.HALT2
            message = (
                f"Missing pairwise statistical column(s): {aggregate_missing_cols}"
            )
        else:
            code = FlagCode.GREEN
            message = f"All pairwise statistical columns exist as follows: {aggregate_target_cols}"

        return {"code": code, "message": message}

    def check_dge_table_statistical_column_constraints(
        self, dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)
        df_dge = pd.read_csv(target_data_asset.path)

        err_msg = ""

        # check all expected statistic columns present
        # pairwise comparison level
        pairwise_comparisons = dataset.metadata.contrasts.columns
        for statistical_prefix, constraints in self.config[
            "pairwise_columns_prefixes"
        ].items():  # type: ignore
            target_cols: list = [
                f"{statistical_prefix}{comparison}"
                for comparison in pairwise_comparisons
            ]
            target_df_subset: pd.DataFrame = df_dge[target_cols]

            # check non null constraint
            if constraints.get("nonNull") and nonNull(target_df_subset) == False:
                err_msg += f"At least one value in columns {target_cols} fails nonNull constraint."
            # check non negative constraint
            if (
                constraints.get("nonNegative")
                and nonNegative(target_df_subset) == False
            ):
                err_msg += f"At least one value in columns {target_cols} fails nonNegative constraint."

        if err_msg:
            code = FlagCode.HALT2
            message = f"Pairwise statistical column(s) failing constraints: {err_msg}"
        else:
            code = FlagCode.GREEN
            message = f"All pairwise statistical columns met constraints: {constraints}"

        return {"code": code, "message": message}

    def check_dge_table_log2fc_within_reason(
        self, dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)
        df_dge = pd.read_csv(target_data_asset.path)

        err_msg = ""
        for comparision in dataset.metadata.contrasts.columns:
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
                    < self.config["log2fc_cross_method_percent_difference_threshold"]
                )
                * 100
            )
            # flag if not enough within tolerance
            if (
                percent_within_tolerance
                < self.config["log2fc_cross_method_tolerance_percent"]
            ):
                err_msg += (
                    f"For comparison: '{comparision}' {percent_within_tolerance:.2f} % of genes have absolute percent differences "
                    f"(between log2fc direct computation and DESeq2's approach) "
                    f"less than {self.config['log2fc_cross_method_percent_difference_threshold']} % which does not met the minimum percentage "
                    f"({self.config['log2fc_cross_method_tolerance_percent']} %) of genes required.  "
                    f"This may indicate misassigned or misaligned columns. "
                )

        if err_msg:
            code = FlagCode.HALT2
            message = f"Pairwise statistical column(s) failing constraints: {err_msg}"
        else:
            code = FlagCode.GREEN
            message = f"All log2fc within reason"

        return {"code": code, "message": message}

    def check_dge_table_group_factor_columns_exist(
        self, dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)

        # read in dataframe
        df_dge = pd.read_csv(target_data_asset.path)

        err_msg = ""

        # factorGroup level
        factorGroups = list(
            set(dataset.metadata.factor_groups.values())
        )  # list-set to dedupe
        for statistical_prefix, constraints in self.config[
            "group_factorwise_columns_prefixes"
        ].items():  # type: ignore
            target_cols = [f"{statistical_prefix}{group}" for group in factorGroups]
            # check existense first and bail if any don't exist
            if missing_cols := set(target_cols) - set(df_dge.columns):
                err_msg += f"Missing groupFactor statistical column(s): {missing_cols}"

        if err_msg:
            code = FlagCode.HALT2
            message = f"Pairwise statistical column(s) failing constraints: {err_msg}"
        else:
            code = FlagCode.GREEN
            message = f"All pairwise statistical columns met constraints: {constraints}"

        return {"code": code, "message": message}

    def check_dge_table_group_factor_column_constraints(
        self, dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):

        err_msg = ""
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)

        # read in dataframe
        df_dge = pd.read_csv(target_data_asset.path)

        # factorGroup level
        factorGroups = list(
            set(dataset.metadata.factor_groups.values())
        )  # list-set to dedupe
        for statistical_prefix, constraints in self.config[
            "group_factorwise_columns_prefixes"
        ].items():  # type: ignore
            target_cols = [f"{statistical_prefix}{group}" for group in factorGroups]
            target_df_subset = df_dge[target_cols]
            # check non null constraint
            if constraints.get("nonNull") and nonNull(target_df_subset) == False:
                err_msg += f"At least one value in columns {target_cols} fails nonNull constraint."
            # check non negative constraint
            if (
                constraints.get("nonNegative")
                and nonNegative(target_df_subset) == False
            ):
                err_msg += f"At least one value in columns {target_cols} fails nonNegative constraint."

        if err_msg:
            code = FlagCode.HALT2
            message = f"Column(s) failing constraints: {err_msg}"
        else:
            code = FlagCode.GREEN
            message = f"All columns met constraints: {constraints}"

        return {"code": code, "message": message}

    def check_dge_table_fixed_statistical_columns_exist(
        self, dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):
        err_msg = ""
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)

        # read in dataframe
        df_dge = pd.read_csv(target_data_asset.path)

        # fixed stat columns level
        for target_col, constraints in self.config["fixed_stats_columns"].items():  # type: ignore
            # check existense first and bail if any don't exist
            if missing_cols := {target_col} - set(df_dge.columns):
                err_msg += f"Missing fixed statistical column(s): {missing_cols}"
                continue
            target_df_subset = df_dge[target_col]

            # check non null constraint
            if constraints.get("nonNull") and nonNull(target_df_subset) == False:
                err_msg += f"At least one value in column ['{target_col}'] fails nonNull constraint."
            # check non negative constraint
            if (
                constraints.get("nonNegative")
                and nonNegative(target_df_subset) == False
            ):
                err_msg += f"At least one value in column ['{target_col}'] fails nonNegative constraint."
        if err_msg:
            code = FlagCode.HALT2
            message = f"Column(s) failing constraints: {err_msg}"
        else:
            code = FlagCode.GREEN
            message = f"All columns met constraints: {constraints}"

        return {"code": code, "message": message}

    def check_dge_table_fixed_statistical_column_constraints(
        self, dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):
        err_msg = ""
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)

        # read in dataframe
        df_dge = pd.read_csv(target_data_asset.path)

        # fixed stat columns level
        for target_col, constraints in self.config["fixed_stats_columns"].items():  # type: ignore
            target_df_subset = df_dge[target_col]

            # check non null constraint
            if constraints.get("nonNull") and nonNull(target_df_subset) == False:
                err_msg += f"At least one value in column ['{target_col}'] fails nonNull constraint."
            # check non negative constraint
            if (
                constraints.get("nonNegative")
                and nonNegative(target_df_subset) == False
            ):
                err_msg += f"At least one value in column ['{target_col}'] fails nonNegative constraint."
        if err_msg:
            code = FlagCode.HALT2
            message = f"Column(s) failing constraints: {err_msg}"
        else:
            code = FlagCode.GREEN
            message = f"All columns met constraints: {constraints}"

        return {"code": code, "message": message}

    def check_dge_table_group_means_computation(
        self, dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):
        err_msg = ""
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)

        # read in dataframe
        df_dge = pd.read_csv(target_data_asset.path)

        # mathematical checks
        groups: list[str] = list(
            {group for group in dataset.metadata.factor_groups.values()}
        )
        aggregate_target_cols: list = list()
        # check means are computed correctly
        for query_group in groups:
            query_column = f"Group.Mean_{query_group}"
            aggregate_target_cols.append(query_column)
            group_samples = [
                sample
                for sample, this_group in dataset.metadata.factor_groups.items()
                if this_group == query_group
            ]
            abs_percent_difference = abs(
                (
                    (
                        (
                            df_dge[group_samples].mean(axis="columns")
                            - df_dge[query_column]
                        )
                        / df_dge[query_column]
                    )
                    * 100
                )
            )
            within_tolerance = abs_percent_difference < self.config["float_tolerance"]
            if not within_tolerance.all() == True:
                err_msg += f"Group Mean value in table is out of float tolerance. This means {query_group} has improperly computed values"

        if err_msg:
            code = FlagCode.HALT2
            message = err_msg
        else:
            code = FlagCode.GREEN
            message = f"All following group mean columns are correctly computed: {aggregate_target_cols}"

        return {"code": code, "message": message}

    def check_viz_pca_table_index_and_columns_exist(
        self,
        dataset: TemplateDataset,
        componentTarget: str,
        dataAssetTarget: str = "visualizationPCATableCSV",
    ) -> str:
        err_msg = ""
        target_component = getattr(dataset, componentTarget)
        target_asset = getattr(target_component, dataAssetTarget)

        # read into dataframe
        df = pd.read_csv(target_asset.path, index_col=0)

        # check all samples included
        if missing_samples := set(dataset.samples.keys()) - set(df.index):
            err_msg += f"Missing samples in index: {missing_samples}"

        # check all expected columns exist
        if missing_cols := set(self.config["expected_vis_pca_columns"]) - set(df.columns):  # type: ignore
            err_msg += f"Missing expected columns: {missing_cols}"

        if err_msg:
            code = FlagCode.HALT2
            message = err_msg
        else:
            code = FlagCode.GREEN
            message = f"PCA Table has all the samples in the index and these columns exist: {list(self.config['expected_vis_pca_columns'])}"

        return {"code": code, "message": message}

    def _viz_output_table_check(
        self, dataset: TemplateDataset, componentTarget: str
    ) -> str:
        """Since this effectively extends the differential expression table,
        run that first and build on the error message as needed"""
        err_msg = self._differential_expression_table_check(
            dataset, componentTarget, componentDataAsset="visualizationTableCSV"
        )

        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, "visualizationTableCSV")

        # read in dataframe
        df = pd.read_csv(target_data_asset.path)

        # check all expected columns exists (all unique to the viz table)
        # check all expected statistic columns present
        # pairwise comparison level
        pairwise_comparisons = dataset.metadata.contrasts.columns
        for statistical_prefix, constraints in self.config[
            "viz_pairwise_columns_prefixes"
        ].items():  # type: ignore
            target_cols: list = [
                f"{statistical_prefix}{comparison}"
                for comparison in pairwise_comparisons
            ]
            # check existense first and bail if any don't exist
            if missing_cols := set(target_cols) - set(df.columns):
                err_msg += f"Missing pairwise statistical column(s): {missing_cols}"
                continue
            target_df_subset: pd.DataFrame = df[target_cols]

            # check non null constraint
            if constraints.get("nonNull") and nonNull(target_df_subset) == False:
                err_msg += f"At least one value in columns {target_cols} fails nonNull constraint."

            # check non negative constraint
            if (
                constraints.get("nonNegative")
                and nonNegative(target_df_subset) == False
            ):
                err_msg += f"At least one value in columns {target_cols} fails nonNegative constraint."

            # check allowed values constraint
            if (
                constraints.get("allowedValues")
                and onlyAllowedValues(
                    target_df_subset, constraints.get("allowedValues")
                )
                == False
            ):
                err_msg += f"At least one value in columns {target_cols} fails allowedValues constraint (allowed values: {constraints.get('allowedValues')})."

        return err_msg


class DATASET_DIFFERENTIALGENEEXPRESSION_ANNOTATED_TABLE_0001(Check):
    USE_SUBCHECK_MONADS = True

    config = {
        "expected_tables": [
            "differential_expression.csv",
            "visualization_output_table.csv",
            "visualization_PCA_table.csv",
        ],
        # Expected column name, but dependent on dataset organism
        "dge_table_master_annotation_keys": {
            "Arabidopsis thaliana": "TAIR",
            "_DEFAULT": "ENSEMBL",
        },
        "dge_table_expected_annotation_columns": [
            "SYMBOL",
            "GENENAME",
            "REFSEQ",
            "ENTREZID",
            "STRING_id",
            "GOSLIM_IDS",
        ],
        # includes column specific constraints
        # these prefix as follows {prefix}{pairWiseFactorGroupComparison}
        "pairwise_columns_prefixes": {
            "Log2fc_": {"nonNull": True},
            "Stat_": {"nonNull": True},
            # can be removed from analysis before p-value and adj-p-value assessed
            # ref: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-are-some-p-values-set-to-na
            "P.value_": {"nonNegative": True, "nonNull": False},
            "Adj.p.value_": {"nonNegative": True, "nonNull": False},
        },
        "viz_pairwise_columns_prefixes": {
            "Log2_Adj.p.value_": {"nonNull": False},
            "Sig.1_": {"allowedValues": [False, True], "nonNull": False},
            "Sig.05_": {"allowedValues": [False, True], "nonNull": False},
            "Log2_P.value_": {"nonNegative": False, "nonNull": False},
            "Updown_": {"allowedValues": [1, 0, -1], "nonNull": True},
        },
        # these prefix as follows {prefix}{FactorGroup}
        "group_factorwise_columns_prefixes": {
            "Group.Mean_": {"nonNull": True, "nonNegative": True},
            "Group.Stdev_": {"nonNull": True, "nonNegative": True},
        },
        "fixed_stats_columns": {
            "All.mean": {"nonNull": True, "nonNegative": True},
            "All.stdev": {"nonNull": True, "nonNegative": True},
            "LRT.p.value": {"nonNull": False, "nonNegative": True},
        },
        "sample_counts_constraints": {"nonNegative": True},
        "expected_vis_pca_columns": [
            "PC1",
            "PC2",
        ],  # more may be included but these are REQUIRED
        "float_tolerance": 0.0001,  # PERCENT
        # TODO: DISCUSS, these baseline values, should indicate a very heavy left-hand skewed histogram of differences - JDO
        "log2fc_cross_method_percent_difference_threshold": 10,  # PERCENT
        "log2fc_cross_method_tolerance_percent": 50,  # PERCENT
        # PERCENT difference minimum between group means to included
        "log2fc_cross_method_sign_check_group_mean_difference_threshold": 50,
        # PERCENT genes allowed sign inversions between methods for groups that meet
        #   log2fc_cross_method_sign_check_group_mean_difference_threshold minimum
        "log2fc_cross_method_sign_check_tolerance_percent": 0,
        # "middle": MIDDLE.median,
        # "yellow_standard_deviation_threshold": 2,
        # "red_standard_deviation_threshold": 4,
    }
    description = (
        "Check that the differential expression outputs exist (source from the deseq2 script) and  "
        "the following tables: {expected_tables}.  "
        "For studies with ERCC spike-in, performs the same check on analogous tables. "
        "Additional performs the file specific validations: "
        "- contrasts.csv: Includes all the existing comparison groups (based on factor values in the metadata) and is formatted correctly"
        "- differential_expression.csv:  Includes expected annotation columns {dge_table_expected_annotation_columns}, includes a master annotation key "
        "column dependent on the dataset organism as follows: {dge_table_master_annotation_keys} ,"
        "includes sample count columns for all samples, all sample count values are non-negative, "
        "all pairwise comparision columns exist with the following prefixes and adhere to the following constraints: {pairwise_columns_prefixes} "
        "all groupFactorWise statistics columns exists with the following prefixes and adhere to the following constraints: {group_factorwise_columns_prefixes} "
        "all fixed statistics columns exist and adhere to the following constraints: {fixed_stats_columns} "
        " - visualization_PCA_table.csv: All samples in index and at the following columns exist {expected_vis_pca_columns} "
        " - visualization_output_table.csv: Performs same checks as differential_expression.csv as well as, "
        "ensuring the additional pairwise comparision columns exist with the following prefixes and "
        "adhere to the following constraints: {expected_vis_pca_columns}. "
        "Confirms that gene counts between differential expression table and normalized counts tables are the same. "
        "Confirms that computations match expectations with respect to following operations: (Float tolerance: +/-{float_tolerance} %)"
        "- Group means are correctly computed from normalized counts "
        "- log2FC values (computed with DESeq2's MLE approach) are comparable to direct computation with log2( mean(group1) / mean(group2) ), specifically "
        "checking if at least {log2fc_cross_method_tolerance_percent} % of all genes have absolute percent differences between methods "
        "less than {log2fc_cross_method_percent_difference_threshold} % "
    )
    flag_desc = {}

    def validate_func(self: Check, dataset: TemplateDataset) -> Flag:
        return (
            Flaggable(
                dataset=dataset,
                componentTarget="differentialGeneExpression",
                componentDataAsset="annotatedTableCSV",
            )
            .bind(self.check_dge_table_annotation_columns_exist)
            .bind(self.check_dge_table_sample_columns_exist)
            .bind(self.check_dge_table_sample_column_constraints)
            .bind(self.check_dge_table_statistical_columns_exist)
            .bind(self.check_dge_table_statistical_column_constraints)
            .bind(self.check_dge_table_group_factor_columns_exist)
            .bind(self.check_dge_table_group_factor_column_constraints)
            .bind(self.check_dge_table_fixed_statistical_columns_exist)
            .bind(self.check_dge_table_fixed_statistical_column_constraints)
            .bind(self.check_dge_table_log2fc_within_reason)
            .bind(self.check_dge_table_group_means_computation)
        ).export_flag_data_to_Flag(check=self)
        all_flags.extend(
            (
                Flaggable(
                    dataset=dataset,
                    componentTarget="differentialGeneExpression",
                    componentDataAsset="visualizationTableCSV",
                )
                .bind(self.check_dge_table_annotation_columns_exist)
                .bind(self.check_dge_table_sample_columns_exist)
                .bind(self.check_dge_table_sample_column_constraints)
                .bind(self.check_dge_table_statistical_columns_exist)
                .bind(self.check_dge_table_statistical_column_constraints)
                .bind(self.check_dge_table_group_factor_columns_exist)
                .bind(self.check_dge_table_group_factor_column_constraints)
                .bind(self.check_dge_table_fixed_statistical_columns_exist)
                .bind(self.check_dge_table_fixed_statistical_column_constraints)
                .bind(self.check_dge_table_log2fc_within_reason)
                .bind(self.check_dge_table_group_means_computation)
                .bind(self.check_viz_pca_table_index_and_columns_exist)
            ).flag_data
        )
        if dataset.metadata.has_ercc:
            all_flags.extend(
                (
                    Flaggable(
                        dataset=dataset,
                        componentTarget="differentialGeneExpressionERCC",
                    )
                    .bind(self.check_contrasts_headers)
                    .bind(self.check_contrasts_rows)
                ).flag_data
            )
            all_flags.extend(
                (
                    Flaggable(
                        dataset=dataset,
                        componentTarget="differentialGeneExpressionERCC",
                        componentDataAsset="annotatedTableCSV",
                    )
                    .bind(self.check_dge_table_annotation_columns_exist)
                    .bind(self.check_dge_table_sample_columns_exist)
                    .bind(self.check_dge_table_sample_column_constraints)
                    .bind(self.check_dge_table_statistical_columns_exist)
                    .bind(self.check_dge_table_statistical_column_constraints)
                    .bind(self.check_dge_table_group_factor_columns_exist)
                    .bind(self.check_dge_table_group_factor_column_constraints)
                    .bind(self.check_dge_table_fixed_statistical_columns_exist)
                    .bind(self.check_dge_table_fixed_statistical_column_constraints)
                    .bind(self.check_dge_table_log2fc_within_reason)
                    .bind(self.check_dge_table_group_means_computation)
                ).flag_data
            )
            all_flags.extend(
                (
                    Flaggable(
                        dataset=dataset,
                        componentTarget="differentialGeneExpressionERCC",
                        componentDataAsset="visualizationTableCSV",
                    )
                    .bind(self.check_dge_table_annotation_columns_exist)
                    .bind(self.check_dge_table_sample_columns_exist)
                    .bind(self.check_dge_table_sample_column_constraints)
                    .bind(self.check_dge_table_statistical_columns_exist)
                    .bind(self.check_dge_table_statistical_column_constraints)
                    .bind(self.check_dge_table_group_factor_columns_exist)
                    .bind(self.check_dge_table_group_factor_column_constraints)
                    .bind(self.check_dge_table_fixed_statistical_columns_exist)
                    .bind(self.check_dge_table_fixed_statistical_column_constraints)
                    .bind(self.check_dge_table_log2fc_within_reason)
                    .bind(self.check_dge_table_group_means_computation)
                    .bind(self.check_viz_pca_table_index_and_columns_exist)
                ).flag_data
            )
        # unwrap and return flag_data
        print(1)

    @staticmethod
    def check_contrasts_headers(
        dataset: TemplateDataset, componentTarget: str
    ) -> FlagEntry:
        # extract target Component
        target_component = getattr(dataset, componentTarget)
        # extract dicts for deseq2 contrasts and the metadata formatted one here
        # make sure to read in explicit index column for deseq2
        dict_deseq2: Dict = pd.read_csv(
            target_component.contrastsCSV.path, index_col=0
        ).to_dict(orient="list")
        dict_data_model: Dict = dataset.metadata.contrasts.to_dict(orient="list")

        # check that all headers are present
        deseq2_headers = set(dict_deseq2.keys())
        data_model_headers = set(dict_data_model.keys())
        if deseq2_headers != data_model_headers:
            code = FlagCode.HALT1
            message = (
                f"Contrast headers do not match expecations "
                f"(based on Factor Values in {dataset.metadata.runsheet.path.name}). "
                f"Mismatched headers: (runsheet-unique: {data_model_headers - deseq2_headers}) "
                f"(contrasts-file-unique: {deseq2_headers - data_model_headers})"
            )
        else:
            code = FlagCode.GREEN
            message = f"Contrast headers correctly formed (based on Factor Values in {dataset.metadata.runsheet.path.name})"
        return {"code": code, "message": message}

    @staticmethod
    def check_contrasts_rows(
        dataset: TemplateDataset, componentTarget: str
    ) -> FlagEntry:
        # extract target Component
        target_component = getattr(dataset, componentTarget)
        # extract dicts for deseq2 contrasts and the metadata formatted one here
        # make sure to read in explicit index column for deseq2
        dict_deseq2: Dict = pd.read_csv(
            target_component.contrastsCSV.path, index_col=0
        ).to_dict(orient="list")
        dict_data_model: Dict = dataset.metadata.contrasts.to_dict(orient="list")

        # check contents of each column matches expecatation (group1 and group2 formatted as expected)
        # this also rechecks headers (keys) but that is caught in the prior validation
        if dict_deseq2 != dict_data_model:
            code = FlagCode.HALT1
            message = (
                f"Rows don't match expectations. Deseq2: {dict_deseq2}. "
                f"DataModel (from metadata source): {dict_data_model}"
            )
        else:
            code = FlagCode.GREEN
            message = f"Contrast rows correctly formed (based on Factor Values in {dataset.metadata.runsheet.path.name})"
        return {"code": code, "message": message}

    def check_dge_table_annotation_columns_exist(
        self,
        dataset: TemplateDataset,
        componentTarget: str,
        componentDataAsset: str = "annotatedTableCSV",
    ):
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)

        # read in dataframe
        df_dge = pd.read_csv(target_data_asset.path)

        # check all constant columns exist
        missing_constant_columns: set
        master_key = self.config["dge_table_master_annotation_keys"].get(
            dataset.metadata.organism,
            self.config["dge_table_master_annotation_keys"]["_DEFAULT"],
        )
        log.debug(
            f"Resolved master annotation key for {dataset.metadata.organism} is {master_key}"
        )
        expected_columns: list = self.config["dge_table_expected_annotation_columns"] + [master_key]  # type: ignore
        if missing_constant_columns := set(expected_columns) - set(df_dge.columns):
            code = FlagCode.HALT2
            err_msg += f"Annotation Columns missing: {missing_constant_columns}"
        else:
            code = FlagCode.GREEN
            message = f"All the following expected annotation columns exist: {expected_columns}"

        return {"code": code, "message": message}

    @staticmethod
    def check_dge_table_sample_columns_exist(
        dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)
        df_dge = pd.read_csv(target_data_asset.path)
        # check all sample counts columns exist
        expected_samples = set(dataset.samples.keys())
        if missing_samples := expected_samples - set(df_dge.columns):
            code = FlagCode.HALT2
            message = f"Sample Count Columns missing: {missing_samples}"
        else:
            code = FlagCode.GREEN
            message = f"Sample Count Columns exist as follows: {list(df_dge.columns)}"

        # check that they met constraints
        # all sample column counts are not negative
        if not (df_dge[list(expected_samples)] >= 0).all(axis=None):
            err_msg += (
                f"Sample Count Columns include negative values: {missing_samples}"
            )

        return {"code": code, "message": message}

    @staticmethod
    def check_dge_table_sample_column_constraints(
        dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)
        df_dge = pd.read_csv(target_data_asset.path)

        expected_samples = set(dataset.samples.keys())

        # check that they met constraints
        # all sample column counts are not negative
        if not (df_dge[list(expected_samples)] >= 0).all(axis=None):
            code = FlagCode.HALT2
            message = f"Sample Count Columns include negative values"
        else:
            code = FlagCode.GREEN
            message = f"Sample Count Columns meet constraints: 'Non-negative'"

        return {"code": code, "message": message}

    def check_dge_table_statistical_columns_exist(
        self, dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)
        df_dge = pd.read_csv(target_data_asset.path)

        aggregate_missing_cols: list = list()
        aggregate_target_cols: list = list()

        # check all expected statistic columns present
        # pairwise comparison level
        pairwise_comparisons = dataset.metadata.contrasts.columns
        for statistical_prefix in self.config[
            "pairwise_columns_prefixes"
        ].keys():  # type: ignore
            target_cols: list = [
                f"{statistical_prefix}{comparison}"
                for comparison in pairwise_comparisons
            ]
            aggregate_target_cols += target_cols
            # check existense first and bail if any don't exist
            if missing_cols := set(target_cols) - set(df_dge.columns):
                aggregate_missing_cols += missing_cols

        if aggregate_missing_cols:
            code = FlagCode.HALT2
            message = (
                f"Missing pairwise statistical column(s): {aggregate_missing_cols}"
            )
        else:
            code = FlagCode.GREEN
            message = f"All pairwise statistical columns exist as follows: {aggregate_target_cols}"

        return {"code": code, "message": message}

    def check_dge_table_statistical_column_constraints(
        self, dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)
        df_dge = pd.read_csv(target_data_asset.path)

        err_msg = ""

        # check all expected statistic columns present
        # pairwise comparison level
        pairwise_comparisons = dataset.metadata.contrasts.columns
        for statistical_prefix, constraints in self.config[
            "pairwise_columns_prefixes"
        ].items():  # type: ignore
            target_cols: list = [
                f"{statistical_prefix}{comparison}"
                for comparison in pairwise_comparisons
            ]
            target_df_subset: pd.DataFrame = df_dge[target_cols]

            # check non null constraint
            if constraints.get("nonNull") and nonNull(target_df_subset) == False:
                err_msg += f"At least one value in columns {target_cols} fails nonNull constraint."
            # check non negative constraint
            if (
                constraints.get("nonNegative")
                and nonNegative(target_df_subset) == False
            ):
                err_msg += f"At least one value in columns {target_cols} fails nonNegative constraint."

        if err_msg:
            code = FlagCode.HALT2
            message = f"Pairwise statistical column(s) failing constraints: {err_msg}"
        else:
            code = FlagCode.GREEN
            message = f"All pairwise statistical columns met constraints: {constraints}"

        return {"code": code, "message": message}

    def check_dge_table_log2fc_within_reason(
        self, dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)
        df_dge = pd.read_csv(target_data_asset.path)

        err_msg = ""
        for comparision in dataset.metadata.contrasts.columns:
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
                    < self.config["log2fc_cross_method_percent_difference_threshold"]
                )
                * 100
            )
            # flag if not enough within tolerance
            if (
                percent_within_tolerance
                < self.config["log2fc_cross_method_tolerance_percent"]
            ):
                err_msg += (
                    f"For comparison: '{comparision}' {percent_within_tolerance:.2f} % of genes have absolute percent differences "
                    f"(between log2fc direct computation and DESeq2's approach) "
                    f"less than {self.config['log2fc_cross_method_percent_difference_threshold']} % which does not met the minimum percentage "
                    f"({self.config['log2fc_cross_method_tolerance_percent']} %) of genes required.  "
                    f"This may indicate misassigned or misaligned columns. "
                )

        if err_msg:
            code = FlagCode.HALT2
            message = f"Pairwise statistical column(s) failing constraints: {err_msg}"
        else:
            code = FlagCode.GREEN
            message = f"All log2fc within reason"

        return {"code": code, "message": message}

    def check_dge_table_group_factor_columns_exist(
        self, dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)

        # read in dataframe
        df_dge = pd.read_csv(target_data_asset.path)

        err_msg = ""

        # factorGroup level
        factorGroups = list(
            set(dataset.metadata.factor_groups.values())
        )  # list-set to dedupe
        for statistical_prefix, constraints in self.config[
            "group_factorwise_columns_prefixes"
        ].items():  # type: ignore
            target_cols = [f"{statistical_prefix}{group}" for group in factorGroups]
            # check existense first and bail if any don't exist
            if missing_cols := set(target_cols) - set(df_dge.columns):
                err_msg += f"Missing groupFactor statistical column(s): {missing_cols}"

        if err_msg:
            code = FlagCode.HALT2
            message = f"Pairwise statistical column(s) failing constraints: {err_msg}"
        else:
            code = FlagCode.GREEN
            message = f"All pairwise statistical columns met constraints: {constraints}"

        return {"code": code, "message": message}

    def check_dge_table_group_factor_column_constraints(
        self, dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):

        err_msg = ""
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)

        # read in dataframe
        df_dge = pd.read_csv(target_data_asset.path)

        # factorGroup level
        factorGroups = list(
            set(dataset.metadata.factor_groups.values())
        )  # list-set to dedupe
        for statistical_prefix, constraints in self.config[
            "group_factorwise_columns_prefixes"
        ].items():  # type: ignore
            target_cols = [f"{statistical_prefix}{group}" for group in factorGroups]
            target_df_subset = df_dge[target_cols]
            # check non null constraint
            if constraints.get("nonNull") and nonNull(target_df_subset) == False:
                err_msg += f"At least one value in columns {target_cols} fails nonNull constraint."
            # check non negative constraint
            if (
                constraints.get("nonNegative")
                and nonNegative(target_df_subset) == False
            ):
                err_msg += f"At least one value in columns {target_cols} fails nonNegative constraint."

        if err_msg:
            code = FlagCode.HALT2
            message = f"Column(s) failing constraints: {err_msg}"
        else:
            code = FlagCode.GREEN
            message = f"All columns met constraints: {constraints}"

        return {"code": code, "message": message}

    def check_dge_table_fixed_statistical_columns_exist(
        self, dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):
        err_msg = ""
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)

        # read in dataframe
        df_dge = pd.read_csv(target_data_asset.path)

        # fixed stat columns level
        for target_col, constraints in self.config["fixed_stats_columns"].items():  # type: ignore
            # check existense first and bail if any don't exist
            if missing_cols := {target_col} - set(df_dge.columns):
                err_msg += f"Missing fixed statistical column(s): {missing_cols}"
                continue
            target_df_subset = df_dge[target_col]

            # check non null constraint
            if constraints.get("nonNull") and nonNull(target_df_subset) == False:
                err_msg += f"At least one value in column ['{target_col}'] fails nonNull constraint."
            # check non negative constraint
            if (
                constraints.get("nonNegative")
                and nonNegative(target_df_subset) == False
            ):
                err_msg += f"At least one value in column ['{target_col}'] fails nonNegative constraint."
        if err_msg:
            code = FlagCode.HALT2
            message = f"Column(s) failing constraints: {err_msg}"
        else:
            code = FlagCode.GREEN
            message = f"All columns met constraints: {constraints}"

        return {"code": code, "message": message}

    def check_dge_table_fixed_statistical_column_constraints(
        self, dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):
        err_msg = ""
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)

        # read in dataframe
        df_dge = pd.read_csv(target_data_asset.path)

        # fixed stat columns level
        for target_col, constraints in self.config["fixed_stats_columns"].items():  # type: ignore
            target_df_subset = df_dge[target_col]

            # check non null constraint
            if constraints.get("nonNull") and nonNull(target_df_subset) == False:
                err_msg += f"At least one value in column ['{target_col}'] fails nonNull constraint."
            # check non negative constraint
            if (
                constraints.get("nonNegative")
                and nonNegative(target_df_subset) == False
            ):
                err_msg += f"At least one value in column ['{target_col}'] fails nonNegative constraint."
        if err_msg:
            code = FlagCode.HALT2
            message = f"Column(s) failing constraints: {err_msg}"
        else:
            code = FlagCode.GREEN
            message = f"All columns met constraints: {constraints}"

        return {"code": code, "message": message}

    def check_dge_table_group_means_computation(
        self, dataset: TemplateDataset, componentTarget: str, componentDataAsset: str
    ):
        err_msg = ""
        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, componentDataAsset)

        # read in dataframe
        df_dge = pd.read_csv(target_data_asset.path)

        # mathematical checks
        groups: list[str] = list(
            {group for group in dataset.metadata.factor_groups.values()}
        )
        aggregate_target_cols: list = list()
        # check means are computed correctly
        for query_group in groups:
            query_column = f"Group.Mean_{query_group}"
            aggregate_target_cols.append(query_column)
            group_samples = [
                sample
                for sample, this_group in dataset.metadata.factor_groups.items()
                if this_group == query_group
            ]
            abs_percent_difference = abs(
                (
                    (
                        (
                            df_dge[group_samples].mean(axis="columns")
                            - df_dge[query_column]
                        )
                        / df_dge[query_column]
                    )
                    * 100
                )
            )
            within_tolerance = abs_percent_difference < self.config["float_tolerance"]
            if not within_tolerance.all() == True:
                err_msg += f"Group Mean value in table is out of float tolerance. This means {query_group} has improperly computed values"

        if err_msg:
            code = FlagCode.HALT2
            message = err_msg
        else:
            code = FlagCode.GREEN
            message = f"All following group mean columns are correctly computed: {aggregate_target_cols}"

        return {"code": code, "message": message}

    def check_viz_pca_table_index_and_columns_exist(
        self,
        dataset: TemplateDataset,
        componentTarget: str,
        dataAssetTarget: str = "visualizationPCATableCSV",
    ) -> str:
        err_msg = ""
        target_component = getattr(dataset, componentTarget)
        target_asset = getattr(target_component, dataAssetTarget)

        # read into dataframe
        df = pd.read_csv(target_asset.path, index_col=0)

        # check all samples included
        if missing_samples := set(dataset.samples.keys()) - set(df.index):
            err_msg += f"Missing samples in index: {missing_samples}"

        # check all expected columns exist
        if missing_cols := set(self.config["expected_vis_pca_columns"]) - set(df.columns):  # type: ignore
            err_msg += f"Missing expected columns: {missing_cols}"

        if err_msg:
            code = FlagCode.HALT2
            message = err_msg
        else:
            code = FlagCode.GREEN
            message = f"PCA Table has all the samples in the index and these columns exist: {list(self.config['expected_vis_pca_columns'])}"

        return {"code": code, "message": message}

    def _viz_output_table_check(
        self, dataset: TemplateDataset, componentTarget: str
    ) -> str:
        """Since this effectively extends the differential expression table,
        run that first and build on the error message as needed"""
        err_msg = self._differential_expression_table_check(
            dataset, componentTarget, componentDataAsset="visualizationTableCSV"
        )

        target_component = getattr(dataset, componentTarget)
        target_data_asset = getattr(target_component, "visualizationTableCSV")

        # read in dataframe
        df = pd.read_csv(target_data_asset.path)

        # check all expected columns exists (all unique to the viz table)
        # check all expected statistic columns present
        # pairwise comparison level
        pairwise_comparisons = dataset.metadata.contrasts.columns
        for statistical_prefix, constraints in self.config[
            "viz_pairwise_columns_prefixes"
        ].items():  # type: ignore
            target_cols: list = [
                f"{statistical_prefix}{comparison}"
                for comparison in pairwise_comparisons
            ]
            # check existense first and bail if any don't exist
            if missing_cols := set(target_cols) - set(df.columns):
                err_msg += f"Missing pairwise statistical column(s): {missing_cols}"
                continue
            target_df_subset: pd.DataFrame = df[target_cols]

            # check non null constraint
            if constraints.get("nonNull") and nonNull(target_df_subset) == False:
                err_msg += f"At least one value in columns {target_cols} fails nonNull constraint."

            # check non negative constraint
            if (
                constraints.get("nonNegative")
                and nonNegative(target_df_subset) == False
            ):
                err_msg += f"At least one value in columns {target_cols} fails nonNegative constraint."

            # check allowed values constraint
            if (
                constraints.get("allowedValues")
                and onlyAllowedValues(
                    target_df_subset, constraints.get("allowedValues")
                )
                == False
            ):
                err_msg += f"At least one value in columns {target_cols} fails allowedValues constraint (allowed values: {constraints.get('allowedValues')})."

        return err_msg


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
        code = FlagCode.HALT1
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
        code = FlagCode.HALT1
        message = f"Missing file: {file.name} expected at {str(file)} "

    return {"code": code, "message": message}


def check_fastqgz_file_contents(file: Path, count_lines_to_check: int) -> FlagEntry:
    # check fastq.gz file looks correct
    if count_lines_to_check == -1:
        count_lines_to_check = float("inf")

    lines_with_issues: list[int] = list()

    # check logic
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
                # update every 2,000,000 reads
                if i % 2_000_000 == 0:
                    log.debug(f"Checked {i} lines for {file}")
                    pass

        if not len(lines_with_issues) == 0:
            code = FlagCode.HALT1
            message = (
                f"Following decompressed fastqGZ lines have issues: {lines_with_issues}"
            )
        else:
            code = FlagCode.GREEN
            message = f"First {count_lines_to_check} lines checked found no issues."
    except (EOFError, gzip.BadGzipFile):
        code = FlagCode.HALT2
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
    output = output.stdout.decode()
    if output == "":
        code = FlagCode.GREEN
        message = f"Samtools quickcheck raised no issues"
    else:
        code = FlagCode.HALT1
        message = f"Samtools quickcheck failed on this file with output: {output}"
    return {"code": code, "message": message}


def check_thresholds(
    component: BaseComponent, mqc_key: str, stat_string: str, thresholds: list[dict]
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
                    code = threshold["code"] if code < threshold["code"] else code

    if code == FlagCode.GREEN:
        message = f"Value: ({value}) did not breech any configured thresholds"
    else:
        message = f"Value: ({value}) breeched configured thresholds"
    return {"code": code, "message": message}


def check_metadata_attributes_exist(
    dataset: TemplateDataset, expected_attrs: list[str]
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
        code = FlagCode.HALT1
        message = f"Missing dataset metadata (source from Runsheet): {missing_metadata_fields}"
    return {"code": code, "message": message}


def check_for_outliers(
    dataset: TemplateDataset,
    sample_component: BaseComponent,
    mqc_module: str,
    mqc_plot: str,
    mqc_keys: list[str],
    thresholds: list[dict],
) -> FlagEntry:
    # assume code is GREEN until outliers detected
    code = FlagCode.GREEN
    # dataframe extraction
    df = dataset.getMQCDataFrame(
        sample_component=sample_component, mqc_module=mqc_module, mqc_plot=mqc_plot
    )

    # track for outliers
    outliers: dict[str, dict[str, dict[str, str]]] = defaultdict(lambda: defaultdict(dict))
    for mqc_key in mqc_keys:
        for threshold in thresholds:
            if threshold['middle_fcn'] == "mean":
                middle = df[mqc_key].mean() 
            elif threshold['middle_fcn'] == "median":
                middle = df[mqc_key].median()
            else:
                raise ValueError(f"Cannot compute middle from supplied middle_fcn name: {threshold['middle_fcn']}. Must supply either 'median' or 'mean'")
            
            # bail if standard deviation == 0
            # e.g. if all values are identical (and thus has no outliers)
            if df[mqc_key].std() == 0:
                continue

            # compute difference
            df_diffs = df[mqc_key]-middle

            # compute as number of standard deviations
            df_diffs_in_std = df_diffs / df[mqc_key].std()

            # add to outlier tracker if over the threshold
            for key, value in df_diffs_in_std.iteritems():
                # if an outlier
                if abs(value) > threshold['stdev_threshold']:
                    # track it
                    outliers[key][mqc_module][mqc_key] = value
                    # elevate code if current code is lower severity
                    if code < threshold['code']:
                        code = threshold['code']
    # check logic
    if code == FlagCode.GREEN:
        message = f"No outliers found for {mqc_keys} in {mqc_plot} part of {mqc_module} multiQC module"
    else:
        message = f"Outliers found in {mqc_module} multiQC module as follows: {outliers}"
    return {"code": code, "message": message, "outliers": outliers}

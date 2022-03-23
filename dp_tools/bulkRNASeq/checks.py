from collections import defaultdict
import copy
import enum
import gzip
import logging
from pathlib import Path
from statistics import mean, median, stdev
from typing import Callable, DefaultDict, Dict, List, Tuple
from dp_tools.components.components import RawReadsComponent
from dp_tools.core.entity_model import (
    DataDir,
    DataFile,
    TemplateDataset,
    TemplateSample,
)

log = logging.getLogger(__name__)

from dp_tools.core.check_model import Check, Flag, FlagCode

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


class SAMPLE_RAWREADS_0001(Check):
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
                    "total_sequences"
                ]
                == sample.rawReverseReads.mqcData["FastQC"]["General_Stats"][
                    "total_sequences"
                ]
            ):
                code = FlagCode.HALT2

        return Flag(
            check=self,
            codes=code,
            message_args={
                "missing_components": missing_components,
                "forward_read_count": sample.rawForwardReads.mqcData["FastQC"][
                    "General_Stats"
                ]["total_sequences"]
                if code == FlagCode.HALT2
                else None,
                "reverse_read_count": sample.rawReverseReads.mqcData["FastQC"][
                    "General_Stats"
                ]["total_sequences"]
                if code == FlagCode.HALT2
                else None,
            },
        )


class SAMPLE_TRIMREADS_0001(SAMPLE_RAWREADS_0001):
    ...


class COMPONENT_RAWREADS_0001(Check):
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
                codes=code,
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
            codes=code,
            message_args={
                "lines_with_issues": lines_with_issues,
                "last_line_checked": i,
                "missing_files": missing_files,
            },
        )


class COMPONENT_TRIMREADS_0001(COMPONENT_RAWREADS_0001):
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


class DATASET_METADATA_0001(Check):
    config = {"expected_metadata_attrs": ["paired_end", "has_ercc"]}
    description = "Checks and reports expected metdata required for processing"
    flag_desc = {
        FlagCode.GREEN: "All expected metadata is accessible and populated. {actual_metadata_fields}",
        FlagCode.HALT1: "Missing expected metadata fields: {missing_metadata_fields}",
    }

    def validate_func(self, dataset: TemplateDataset) -> Flag:
        # assume green unless flag condition met
        code = FlagCode.GREEN

        # set up tracker for expected attributes values
        tracked_metadata = dict()
        # and a tracker for missing attributes
        missing_metadata_fields = list()

        for attr in self.config["expected_metadata_attrs"]:
            attr_value = getattr(dataset.metadata, attr, None)
            if attr_value != None:
                tracked_metadata[attr] = attr_value
            else:
                missing_metadata_fields.append(attr)

        # check if any missing_metadata_fields are present
        if missing_metadata_fields:
            code = FlagCode.HALT1

        return Flag(
            check=self,
            codes=code,
            message_args={
                "actual_metadata_fields": tracked_metadata,
                "missing_metadata_fields": missing_metadata_fields,
            },
        )


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
        "yellow_minimum_dominant_strandedness": 75,  # percents
        "halt_minimum_dominant_strandedness": 65,  # percents
    }
    description = (
        "Check that the rseqc analysis stats (source from the rseqc logs) have no outlier values among samples "
        "for the following plots: {plots_all} (Paired end only: {plot_paired_end}). "
        "Yellow Flagged Outliers are defined as a being {yellow_standard_deviation_threshold} - {red_standard_deviation_threshold} standard "
        "deviations away from the {middle.name}. "
        "Red Flagged Outliers are defined as a being {red_standard_deviation_threshold}+ standard "
        "deviations away from the {middle.name}. "
        "Additionally the following is assessed for infer experiment strandedess metrics: "
        "A Yellow Flag is raised in the case that the dominant strandessness is between {yellow_minimum_dominant_strandedness} - {halt_minimum_dominant_strandedness} "
        "A Halt Flag is raised in the case that the dominant strandessness is below {halt_minimum_dominant_strandedness} "
        "Note: the 'dominant strandedness' is the max(datasetwide_average(antisense), datasetwide_average(sense)) "
    )
    flag_desc = {
        FlagCode.GREEN: "No rseqc analysis metric outliers detected for {metrics}",
        FlagCode.YELLOW1: "Outliers detected as follows (values are rounded number of standard deviations from middle): {formatted_outliers}",
        FlagCode.YELLOW2: "The dominant strandedness is {dominant_strandedness}, this is lower than the yellow flag threshold.",
        FlagCode.RED1: "Outliers detected as follows (values are rounded number of standard deviations from middle): {formatted_outliers}",
        FlagCode.HALT1: "The dominant strandedness is {dominant_strandedness}, this is lower than the halting flag threshold.",
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

        # dominant strandedness related subcheck
        # assess dominant strandedness
        yellow_minimum_dominant_strandedness = self.config[
            "yellow_minimum_dominant_strandedness"
        ]
        halt_minimum_dominant_strandedness = self.config[
            "halt_minimum_dominant_strandedness"
        ]

        def get_dominant_strandedness(dataset: TemplateDataset) -> tuple[str, float]:
            df = dataset.getMQCDataFrame(
                sample_component="rSeQCAnalysis",
                mqc_module="RSeQC",
                mqc_plot="Infer experiment",
            )
            dominant_strandedness_name = df.mean()[
                ["Antisense (% Tags)", "Sense (% Tags)"]
            ].idxmax()
            dominant_strandedness_value = df.mean()[
                ["Antisense (% Tags)", "Sense (% Tags)"]
            ].max()
            return (dominant_strandedness_name, dominant_strandedness_value)

        dominant_strandedness = get_dominant_strandedness(dataset)

        # flag based on thresholds
        dominant_strandedness_value = dominant_strandedness[1]
        if (
            yellow_minimum_dominant_strandedness
            > dominant_strandedness_value
            > halt_minimum_dominant_strandedness
        ):
            codes.add(FlagCode.YELLOW2)
        elif dominant_strandedness_value < halt_minimum_dominant_strandedness:
            codes.add(FlagCode.HALT1)

        # set code to green if no flags added

        return Flag(
            codes=codes,
            check=self,
            message_args={
                "outliers": outliers,
                "formatted_outliers": pformat(outliers, formatfloat),
                "dominant_strandedness": dominant_strandedness,
            },
        )

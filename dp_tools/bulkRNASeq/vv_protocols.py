from collections import defaultdict
import enum
import os
from pathlib import Path
from typing import Dict, List, Set, Tuple, Union
import logging

from dp_tools.components.components import (
    DatasetGeneCounts,
    DifferentialGeneExpression,
    GeneCounts,
    GenomeAlignments,
    NormalizedGeneCounts,
    RSeQCAnalysis,
)

log = logging.getLogger(__name__)

from dp_tools.bulkRNASeq.entity import BulkRNASeqDataset, BulkRNASeqSample
from dp_tools.core.check_model import Flag, VVProtocol, ValidationProtocol, FlagCode
from dp_tools.bulkRNASeq.checks import (
    check_fastqgz_file_contents,
    check_forward_and_reverse_reads_counts_match,
    check_file_exists,
    check_bam_file_integrity,
    check_thresholds,
    check_metadata_attributes_exist,
    check_for_outliers,
    MIDDLE
)
from dp_tools.core.entity_model import TemplateComponent
from dp_tools.components import (
    BulkRNASeqMetadataComponent,
    RawReadsComponent,
    TrimReadsComponent,
)


class STAGE(enum.Enum):
    Demultiplexed = 0
    Reads_PreProcessed = 1
    GenomeAligned = 2
    RSeQCAnalysis = 2.01
    GeneCounted = 3
    DGE = 4

    # allow comparing stages
    # used to check if a sufficient stage has been achieved
    def __ge__(self, other):
        return self.value >= other.value

    def __le__(self, other):
        return self.value <= other.value

    @classmethod
    def get_all_preceeding(cls, query_stage):
        return {stage for stage in cls if stage <= query_stage}


def validate_bulkRNASeq(
    dataset: BulkRNASeqDataset, report_args: dict = None, protocol_args: dict = None
):
    if report_args is None:
        report_args = dict()
    if protocol_args is None:
        protocol_args = dict()
    # init validation protocol
    vp = ValidationProtocol(**protocol_args)

    with vp.component_start(
        name=dataset.name,
        description="Validate processing from trim reads through differential gene expression output",
    ):
        # fmt: off
        for sample in dataset.samples.values():
            with vp.component_start(name=sample.name, description="Samples level checks"):
                with vp.component_start(name="Raw Reads", description="Raw reads"):
                    with vp.payload(
                        payloads=([
                            {"file": getattr(getattr(sample.trimForwardReads, "fastqGZ", None), "path", None)},
                            {"file": getattr(getattr(sample.trimReverseReads, "fastqGZ", None), "path", None)},
                        ] if dataset.metadata.paired_end else [
                            {"file": getattr(getattr(sample.trimReads, "fastqGZ", None), "path", None)}
                        ]
                        )
                    ) as RUN:
                        RUN(check_file_exists, description="Check reads files exist")
                        RUN(check_fastqgz_file_contents, config={"count_lines_to_check":200_000_000})

                    with vp.payload(
                        payloads=[
                            {
                                "fwd_reads": sample.trimForwardReads,
                                "rev_reads": sample.trimReverseReads,
                            },
                        ],
                        skip=(not dataset.metadata.paired_end),
                    ) as RUN:
                        RUN(check_forward_and_reverse_reads_counts_match)

                    with vp.component_start(name="multiQC", description="MultiQC checks"):
                        with vp.payload(
                            payloads=([
                                    {"file": getattr(getattr(sample.trimForwardReads, "fastQCmultiQCDirZIP", None), "path", None)},
                                    {"file": getattr(getattr(sample.trimForwardReads, "fastqcReportHTML", None), "path", None)},
                                    {"file": getattr(getattr(sample.trimForwardReads, "fastqcReportZIP", None), "path", None)},
                                    {"file": getattr(getattr(sample.trimReverseReads, "fastQCmultiQCDirZIP", None), "path", None)},
                                    {"file": getattr(getattr(sample.trimReverseReads, "fastqcReportHTML", None), "path", None)},
                                    {"file": getattr(getattr(sample.trimReverseReads, "fastqcReportZIP", None), "path", None)},
                            ] if dataset.metadata.paired_end else [
                                    {"file": getattr(getattr(sample.trimReads, "fastQCmultiQCDirZIP", None), "path", None)},
                                    {"file": getattr(getattr(sample.trimReads, "fastqcReportHTML", None), "path", None)},
                                    {"file": getattr(getattr(sample.trimReads, "fastqcReportZIP", None), "path", None)},
                            ]
                            )
                        ) as RUN:
                            RUN(check_file_exists)

                with vp.component_start(name="Trimmed Reads", description="Trimmed reads"):
                    with vp.payload(
                        payloads=([
                            {"file": getattr(getattr(sample.trimForwardReads, "fastqGZ", None), "path", None)},
                            {"file": getattr(getattr(sample.trimReverseReads, "fastqGZ", None), "path", None)},
                        ] if dataset.metadata.paired_end else [
                            {"file": getattr(getattr(sample.trimReads, "fastqGZ", None), "path", None)}
                        ]
                        )
                    ) as RUN:
                        RUN(check_file_exists, description="Check reads files exist")
                        RUN(check_fastqgz_file_contents, config={"count_lines_to_check":200_000_000})

                    with vp.payload(
                        payloads=[
                            {
                                "fwd_reads": sample.trimForwardReads,
                                "rev_reads": sample.trimReverseReads,
                            },
                        ],
                        skip=(not dataset.metadata.paired_end),
                    ) as RUN:
                        RUN(check_forward_and_reverse_reads_counts_match)

                    with vp.component_start(name="multiQC", description="MultiQC checks"):
                        with vp.payload(
                            payloads=([
                                    {"file": getattr(getattr(sample.trimForwardReads, "fastQCmultiQCDirZIP", None), "path", None)},
                                    {"file": getattr(getattr(sample.trimForwardReads, "fastqcReportHTML", None), "path", None)},
                                    {"file": getattr(getattr(sample.trimForwardReads, "fastqcReportZIP", None), "path", None)},
                                    {"file": getattr(getattr(sample.trimReverseReads, "fastQCmultiQCDirZIP", None), "path", None)},
                                    {"file": getattr(getattr(sample.trimReverseReads, "fastqcReportHTML", None), "path", None)},
                                    {"file": getattr(getattr(sample.trimReverseReads, "fastqcReportZIP", None), "path", None)},
                            ] if dataset.metadata.paired_end else [
                                    {"file": getattr(getattr(sample.trimReads, "fastQCmultiQCDirZIP", None), "path", None)},
                                    {"file": getattr(getattr(sample.trimReads, "fastqcReportHTML", None), "path", None)},
                                    {"file": getattr(getattr(sample.trimReads, "fastqcReportZIP", None), "path", None)},
                            ]
                            )
                        ) as RUN:
                            RUN(check_file_exists)

                    with vp.component_start(
                        name="Trimming Reports", description="Trimming Reports as output by Trim Galore!"
                    ):
                        with vp.payload(
                            payloads=([
                                    {"file": getattr(getattr(sample.trimForwardReads, "trimmingReportTXT", None), "path", None)},
                                    {"file": getattr(getattr(sample.trimReverseReads, "trimmingReportTXT", None), "path", None)},
                            ] if dataset.metadata.paired_end else [
                                    {"file": getattr(getattr(sample.trimReads, "trimmingReportTXT", None), "path", None)},
                            ]
                            )
                        ) as RUN:
                            RUN(check_file_exists, description="Check that Trim Galore reports exist")

                with vp.component_start(name="STAR Alignments", description="STAR Alignment outputs"):
                    with vp.payload(payloads=[
                        {"file": sample.genomeAlignments.alignedToTranscriptomeBam.path},
                        {"file": sample.genomeAlignments.alignedSortedByCoordBam.path},
                        {"file": sample.genomeAlignments.logFinal.path},
                        {"file": sample.genomeAlignments.logProgress.path},
                        {"file": sample.genomeAlignments.logFull.path},
                        {"file": sample.genomeAlignments.sjTab.path},
                    ]) as RUN:
                        RUN(check_file_exists, description="Check that all expected output from STAR is generated")
                    
                    with vp.payload(payloads=[
                        {"file": sample.genomeAlignments.alignedToTranscriptomeBam.path},
                        {"file": sample.genomeAlignments.alignedSortedByCoordBam.path},
                    ]) as RUN:
                        #RUN(check_bam_file_integrity, config={"samtools_bin":#TODO: fill with private config file})
                        RUN(check_bam_file_integrity, config={"samtools_bin":os.environ["SAMTOOLS_BIN"]})

                    with vp.payload(payloads=[
                        {"component": sample.genomeAlignments, "mqc_key":"STAR"},
                    ]) as RUN:
                        # fmt: on
                        RUN(check_thresholds, 
                            config={
                                "stat_string":"uniquely_mapped_percent + multimapped_percent",
                                "thresholds":[
                                            {"code": FlagCode.YELLOW1, "type": "lower", "value": 70},
                                            {"code": FlagCode.RED1, "type": "lower", "value": 50},
                                        ],
                                },
                            description="Check that mapping rates are reasonable, specifically most reads map to the target genome"
                        )
                        RUN(check_thresholds, 
                            config={
                                "stat_string":"multimapped_toomany_percent + multimapped_percent",
                                "thresholds":[
                                            {"code": FlagCode.YELLOW1, "type": "lower", "value": 30},
                                            {"code": FlagCode.RED1, "type": "lower", "value": 15},
                                        ],
                                },
                            description="Check that mapping rates are reasonable, specifically that a considerable amount of reads multimap to the target genome"
                        )

        with vp.component_start(name="Metadata", description="Metadata file validation"):
            with vp.payload(payloads=[
                {"dataset":dataset}
            ]) as RUN:
                RUN(check_metadata_attributes_exist, config={"expected_attrs":["paired_end","has_ercc"]})

        with vp.component_start(name="Raw Reads", description="Raw Reads Outliers Detection"):
            with vp.payload(payloads=[
                {"dataset":dataset, 
                "sample_component":"rawReads",
                "mqc_module":"FastQC"}
            ] if not dataset.metadata.paired_end else [
                {"dataset":dataset, 
                "sample_component":"rawForwardReads",
                "mqc_module":"FastQC"},
                {"dataset":dataset, 
                "sample_component":"rawReverseReads",
                "mqc_module":"FastQC"},
            ]) as RUN:
                RUN(check_for_outliers,
                    config={
                        "mqc_plot":"general_stats",
                        "mqc_keys":[
                            "percent_gc",
                            "avg_sequence_length",
                            "total_sequences",
                            "percent_duplicates",
                        ],
                        "thresholds":[
                            {"code": FlagCode.YELLOW1, "stdev_threshold":0.2, "middle_fcn":'median'},
                            {"code": FlagCode.RED1, "stdev_threshold":0.6, "middle_fcn":'median'},
                        ]
                    }
                )

        with vp.component_start(name="Trim Reads", description="Trimmed Reads Outliers Detection"):
            with vp.payload(payloads=[
                {"dataset":dataset, 
                "sample_component":"trimReads",
                "mqc_module":"FastQC"}
            ] if not dataset.metadata.paired_end else [
                {"dataset":dataset, 
                "sample_component":"trimForwardReads",
                "mqc_module":"FastQC"},
                {"dataset":dataset, 
                "sample_component":"trimReverseReads",
                "mqc_module":"FastQC"},
            ]) as RUN:
                RUN(check_for_outliers,
                    config={
                        "mqc_plot":"general_stats",
                        "mqc_keys":[
                            "percent_gc",
                            "avg_sequence_length",
                            "total_sequences",
                            "percent_duplicates",
                        ],
                        "thresholds":[
                            {"code": FlagCode.YELLOW1, "stdev_threshold":0.2, "middle_fcn":'median'},
                            {"code": FlagCode.RED1, "stdev_threshold":0.6, "middle_fcn":'median'},
                        ]
                    }
                )


    # return report
    return vp.report(**report_args)

'''
class BulkRNASeq_VVProtocol(VVProtocol):
    expected_dataset_class = BulkRNASeqDataset

    # init all checks
    STAGES = STAGE

    # TODO: Move generalized functionality to abc init
    def __init__(
        self,
        dataset: BulkRNASeqDataset,
        config: Union[Tuple[str, str], Path] = ("bulkRNASeq", "Latest"),
        **kwargs,
    ):
        super().__init__(dataset=dataset, config=config, **kwargs)

    def validation_procedure(self):
        samples_payloads = [
            {"sample": sample} for sample in self.dataset.samples.values()
        ]
        if STAGE.Demultiplexed in self._stages_loaded:
            self.batch_check(
                level="dataset",
                payloads=[{"dataset": self.dataset}],
                run=[DATASET_METADATA_0001, DATASET_RAWREADS_0001],
            )
            self.batch_check(
                level="sample", payloads=samples_payloads, run=[SAMPLE_RAWREADS_0001]
            )
            self.batch_check(
                level="component", payloads=samples_payloads, run=[SAMPLE_RAWREADS_0001]
            )
            flags[component].append(
                self.run_check(COMPONENT_RAWREADS_0001, component=component)
            )
        if STAGE.Reads_PreProcessed in self._stages_loaded:
            self.batch_check(
                level="dataset",
                payloads=[{"dataset": self.dataset}],
                run=[DATASET_TRIMREADS_0001],
            )
            self.batch_check(
                level="sample", payloads=samples_payloads, run=[SAMPLE_TRIMREADS_0001]
            )
        if STAGE.GenomeAligned in self._stages_loaded:
            self.batch_check(
                level="dataset",
                payloads=[{"dataset": self.dataset}],
                run=[DATASET_GENOMEALIGNMENTS_0001],
            )
        if STAGE.RSeQCAnalysis in self._stages_loaded:
            self.batch_check(
                level="dataset",
                payloads=[{"dataset": self.dataset}],
                run=[DATASET_RSEQCANALYSIS_0001],
            )
        if STAGE.GeneCounted in self._stages_loaded:
            self.batch_check(
                level="dataset",
                payloads=[{"dataset": self.dataset}],
                run=[DATASET_GENECOUNTS_0001],
            )
        if STAGE.DGE in self._stages_loaded:
            self.batch_check(
                level="dataset",
                payloads=[{"dataset": self.dataset}],
                run=[
                    DATASET_DIFFERENTIALGENEEXPRESSION_CONTRASTS_0001,
                    DATASET_DIFFERENTIALGENEEXPRESSION_ANNOTATED_TABLE_0001,
                ],
            )

    def validate_components(self) -> Dict[TemplateComponent, List[Flag]]:
        flags: Dict[TemplateComponent, List[Flag]] = defaultdict(list)
        # iterate through all components by level
        for component in self.dataset.all_non_empty_components_recursive:
            log.debug(f"Validating component: {component.__class__.__name__}")
            match component:
                case RawReadsComponent():
                    if STAGE.Demultiplexed in self._stages_loaded:
                        flags[component].append(
                            self.run_check(COMPONENT_RAWREADS_0001, component=component)
                        )
                case TrimReadsComponent():
                    if STAGE.Reads_PreProcessed in self._stages_loaded:
                        flags[component].append(
                            self.run_check(
                                COMPONENT_TRIMREADS_0001, component=component
                            )
                        )
                case BulkRNASeqMetadataComponent():
                    flags[component] = list()
                case GenomeAlignments():
                    if STAGE.GenomeAligned in self._stages_loaded:
                        flags[component].append(
                            self.run_check(
                                COMPONENT_GENOMEALIGNMENTS_0001,
                                component=component,
                            )
                        )
                case RSeQCAnalysis():
                    pass  # no component level checks implemented
                case GeneCounts():
                    pass  # no component level checks implemented
                case DatasetGeneCounts():
                    pass  # no component level checks implemented
                case NormalizedGeneCounts():
                    pass  # no component level checks implemented
                case DifferentialGeneExpression():
                    pass  # no component level checks implemented
                case _:
                    raise TypeError(
                        f"Encountered unhandled component type in VV: {component} with type {type(component)}"
                    )
        return flags
'''
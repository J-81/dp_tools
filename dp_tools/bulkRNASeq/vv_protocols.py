from collections import defaultdict
import enum
import os
from pathlib import Path
from typing import Dict, List, Set, Tuple, Union
import yaml
import logging

log = logging.getLogger(__name__)

from dp_tools.bulkRNASeq.entity import BulkRNASeqDataset, BulkRNASeqSample
from dp_tools.core.check_model import ValidationProtocol, FlagCode
from dp_tools.bulkRNASeq.checks import (
    check_fastqgz_file_contents,
    check_forward_and_reverse_reads_counts_match,
    check_file_exists,
    check_bam_file_integrity,
    check_thresholds,
    check_metadata_attributes_exist,
    check_for_outliers,
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
    dataset: BulkRNASeqDataset,
    config_path: Path,
    report_args: dict = None,
    protocol_args: dict = None,
    defer_run: bool = False,
) -> Union[ValidationProtocol, ValidationProtocol.Report]:

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

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

        with vp.component_start(
            name="Metadata", description="Metadata file validation"
        ):
            with vp.payload(payloads=[{"dataset": dataset}]) as ADD:
                ADD(
                    check_metadata_attributes_exist,
                    config=config["Metadata-check_metadata_attributes_exist"],
                )

        with vp.component_start(
            name="Raw Reads", description="Raw Reads Outliers Detection"
        ):
            with vp.payload(
                payloads=[{"dataset": dataset, "sample_component": "rawReads"}]
                if not dataset.metadata.paired_end
                else [
                    {"dataset": dataset, "sample_component": "rawForwardReads"},
                    {"dataset": dataset, "sample_component": "rawReverseReads"},
                ]
            ) as ADD:
                ADD(check_for_outliers, config=config["Raw Reads-check_for_outliers"])

        with vp.component_start(
            name="Trim Reads", description="Trimmed Reads Outliers Detection"
        ):
            with vp.payload(
                payloads=[
                    {
                        "dataset": dataset,
                        "sample_component": "trimReads",
                        "mqc_module": "FastQC",
                    }
                ]
                if not dataset.metadata.paired_end
                else [
                    {
                        "dataset": dataset,
                        "sample_component": "trimForwardReads",
                        "mqc_module": "FastQC",
                    },
                    {
                        "dataset": dataset,
                        "sample_component": "trimReverseReads",
                        "mqc_module": "FastQC",
                    },
                ]
            ) as ADD:
                ADD(
                    check_for_outliers,
                    config={
                        "mqc_plot": "general_stats",
                        "mqc_keys": [
                            "percent_gc",
                            "avg_sequence_length",
                            "total_sequences",
                            "percent_duplicates",
                        ],
                        "thresholds": [
                            {
                                "code": FlagCode.YELLOW1,
                                "stdev_threshold": 0.2,
                                "middle_fcn": "median",
                            },
                            {
                                "code": FlagCode.RED1,
                                "stdev_threshold": 0.6,
                                "middle_fcn": "median",
                            },
                        ],
                    },
                )

        for sample in dataset.samples.values():
            with vp.component_start(
                name=sample.name, description="Samples level checks"
            ):
                with vp.component_start(
                    name="Raw Reads By Sample", description="Raw reads"
                ):
                    with vp.payload(
                        payloads=(
                            [
                                {"file": lambda: sample.rawForwardReads.fastqGZ.path},
                                {"file": lambda: sample.rawReverseReads.fastqGZ.path},
                            ]
                            if dataset.metadata.paired_end
                            else [{"file": lambda: sample.rawReads.fastqGZ.path}]
                        )
                    ) as ADD:
                        ADD(check_file_exists, description="Check reads files exist")
                        ADD(
                            check_fastqgz_file_contents,
                            config={"count_lines_to_check": 200_000_000},
                        )

                    with vp.payload(
                        payloads=[
                            {
                                "fwd_reads": lambda: sample.rawForwardReads,
                                "rev_reads": lambda: sample.rawReverseReads,
                            },
                        ],
                    ) as ADD:
                        ADD(
                            check_forward_and_reverse_reads_counts_match,
                            skip=(not dataset.metadata.paired_end),
                        )

                    with vp.component_start(
                        name="multiQC", description="MultiQC checks"
                    ):
                        with vp.payload(
                            payloads=(
                                [
                                    {
                                        "file": lambda: sample.rawForwardReads.fastQCmultiQCDirZIP.path
                                    },
                                    {
                                        "file": lambda: sample.rawForwardReads.fastqcReportHTML.path
                                    },
                                    {
                                        "file": lambda: sample.rawForwardReads.fastqcReportZIP.path
                                    },
                                    {
                                        "file": lambda: sample.rawReverseReads.fastQCmultiQCDirZIP.path
                                    },
                                    {
                                        "file": lambda: sample.rawReverseReads.fastqcReportHTML.path
                                    },
                                    {
                                        "file": lambda: sample.rawReverseReads.fastqcReportZIP.path
                                    },
                                ]
                                if dataset.metadata.paired_end
                                else [
                                    {
                                        "file": lambda: sample.rawReads.fastQCmultiQCDirZIP.path
                                    },
                                    {
                                        "file": lambda: sample.rawReads.fastqcReportHTML.path
                                    },
                                    {
                                        "file": lambda: sample.rawReads.fastqcReportZIP.path
                                    },
                                ]
                            )
                        ) as ADD:
                            ADD(check_file_exists)

                with vp.component_start(
                    name="Trimmed Reads By Sample", description="Trimmed reads"
                ):
                    with vp.payload(
                        payloads=(
                            [
                                {"file": lambda: sample.trimForwardReads.fastqGZ.path},
                                {"file": lambda: sample.trimReverseReads.fastqGZ.path},
                            ]
                            if dataset.metadata.paired_end
                            else [{"file": lambda: sample.trimReads.fastqGZ.path}]
                        )
                    ) as ADD:
                        ADD(check_file_exists, description="Check reads files exist")
                        ADD(
                            check_fastqgz_file_contents,
                            config={"count_lines_to_check": 200_000_000},
                        )

                    with vp.payload(
                        payloads=[
                            {
                                "fwd_reads": sample.trimForwardReads,
                                "rev_reads": sample.trimReverseReads,
                            },
                        ],
                    ) as ADD:
                        ADD(
                            check_forward_and_reverse_reads_counts_match,
                            skip=(not dataset.metadata.paired_end),
                        )

                    with vp.component_start(
                        name="multiQC", description="MultiQC checks"
                    ):
                        with vp.payload(
                            payloads=(
                                [
                                    {
                                        "file": lambda: sample.trimForwardReads.fastQCmultiQCDirZIP.path
                                    },
                                    {
                                        "file": lambda: sample.trimForwardReads.fastqcReportHTML.path
                                    },
                                    {
                                        "file": lambda: sample.trimForwardReads.fastqcReportZIP.path
                                    },
                                    {
                                        "file": lambda: sample.trimReverseReads.fastQCmultiQCDirZIP.path
                                    },
                                    {
                                        "file": lambda: sample.trimReverseReads.fastqcReportHTML.path
                                    },
                                    {
                                        "file": lambda: sample.trimReverseReads.fastqcReportZIP.path
                                    },
                                ]
                                if dataset.metadata.paired_end
                                else [
                                    {
                                        "file": lambda: sample.trimReads.fastQCmultiQCDirZIP.path
                                    },
                                    {
                                        "file": lambda: sample.trimReads.fastqcReportHTML.path
                                    },
                                    {
                                        "file": lambda: sample.trimReads.fastqcReportZIP.path
                                    },
                                ]
                            )
                        ) as ADD:
                            ADD(check_file_exists)

                    with vp.component_start(
                        name="Trimming Reports",
                        description="Trimming Reports as output by Trim Galore!",
                    ):
                        with vp.payload(
                            payloads=(
                                [
                                    {
                                        "file": lambda: sample.trimForwardReads.trimmingReportTXT.path
                                    },
                                    {
                                        "file": lambda: sample.trimReverseReads.trimmingReportTXT.path
                                    },
                                ]
                                if dataset.metadata.paired_end
                                else [
                                    {
                                        "file": lambda: sample.trimReads.trimmingReportTXT.path
                                    },
                                ]
                            )
                        ) as ADD:
                            ADD(
                                check_file_exists,
                                description="Check that Trim Galore reports exist",
                            )

                with vp.component_start(
                    name="STAR Alignments By Sample",
                    description="STAR Alignment outputs",
                ):
                    with vp.payload(
                        payloads=[
                            {
                                "file": lambda: sample.genomeAlignments.alignedToTranscriptomeBam.path
                            },
                            {
                                "file": lambda: sample.genomeAlignments.alignedSortedByCoordBam.path
                            },
                            {"file": lambda: sample.genomeAlignments.logFinal.path},
                            {"file": lambda: sample.genomeAlignments.logProgress.path},
                            {"file": lambda: sample.genomeAlignments.logFull.path},
                            {"file": lambda: sample.genomeAlignments.sjTab.path},
                        ]
                    ) as ADD:
                        ADD(
                            check_file_exists,
                            description="Check that all expected output from STAR is generated",
                        )

                    with vp.payload(
                        payloads=[
                            {
                                "file": lambda: sample.genomeAlignments.alignedToTranscriptomeBam.path
                            },
                            {
                                "file": lambda: sample.genomeAlignments.alignedSortedByCoordBam.path
                            },
                        ]
                    ) as ADD:
                        # ADD(check_bam_file_integrity, config={"samtools_bin":#TODO: fill with private config file})
                        ADD(
                            check_bam_file_integrity,
                            config={"samtools_bin": os.environ["SAMTOOLS_BIN"]},
                        )

                    with vp.payload(
                        payloads=[
                            {
                                "component": lambda: sample.genomeAlignments,
                                "mqc_key": "STAR",
                            },
                        ]
                    ) as ADD:
                        # fmt: on
                        ADD(
                            check_thresholds,
                            config={
                                "stat_string": "uniquely_mapped_percent + multimapped_percent",
                                "thresholds": [
                                    {
                                        "code": FlagCode.YELLOW1,
                                        "type": "lower",
                                        "value": 70,
                                    },
                                    {
                                        "code": FlagCode.RED1,
                                        "type": "lower",
                                        "value": 50,
                                    },
                                ],
                            },
                            description="Check that mapping rates are reasonable, specifically most reads map to the target genome",
                        )
                        ADD(
                            check_thresholds,
                            config={
                                "stat_string": "multimapped_toomany_percent + multimapped_percent",
                                "thresholds": [
                                    {
                                        "code": FlagCode.YELLOW1,
                                        "type": "lower",
                                        "value": 30,
                                    },
                                    {
                                        "code": FlagCode.RED1,
                                        "type": "lower",
                                        "value": 15,
                                    },
                                ],
                            },
                            description="Check that mapping rates are reasonable, specifically that a considerable amount of reads multimap to the target genome",
                        )

    # return protocol object without running or generating a report
    if defer_run:
        return vp

    vp.run()

    # return report
    return vp.report(**report_args)

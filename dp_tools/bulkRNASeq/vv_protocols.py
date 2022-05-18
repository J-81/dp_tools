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
    check_genebody_coverage_output,
    check_inner_distance_output,
    check_strandedness_assessable_from_infer_experiment,
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
    # fmt: off
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
                payloads=[{"dataset": dataset, "sample_component": "trimReads"}]
                if not dataset.metadata.paired_end
                else [
                    {"dataset": dataset, "sample_component": "trimForwardReads"},
                    {"dataset": dataset, "sample_component": "trimReverseReads"},
                ]
            ) as ADD:
                ADD(
                    check_for_outliers,
                    config=config["Trim Reads-check_for_outliers"],
                )

        with vp.component_start(
            name="STAR Alignments",
            description="Dataset wide checks including outliers detection",
            ):
            with vp.payload(payloads=[
                {"dataset":dataset, "sample_component": "genomeAlignments"}
            ]) as ADD:
                ADD(check_for_outliers, config=config["STAR Alignments-check_for_outliers"])

        with vp.component_start(
            name="RSeQC",
            description="RSeQC submodule outliers checking and other submodule specific dataset wide checks",
        ):
            with vp.payload(
                payloads=[{"dataset": dataset, "sample_component": "rSeQCAnalysis"}]
            ) as ADD:
                ADD(check_for_outliers, description="Check for outliers in geneBody Coverage",config=config["RSeQC-check_for_outliers-geneBody_coverage"])
                ADD(check_for_outliers, description="Check for outliers in infer experiment",config=config["RSeQC-check_for_outliers-infer_experiment"])
                ADD(check_for_outliers, description="Check for outliers in inner distance",config=config["RSeQC-check_for_outliers-inner_distance"], skip=(not dataset.metadata.paired_end))
                ADD(check_for_outliers, description="Check for outliers in read distribution",config=config["RSeQC-check_for_outliers-read_distribution"])

            with vp.payload(payloads=[
                    {"dataset": dataset}
                ]) as ADD:
                ADD(check_strandedness_assessable_from_infer_experiment, config=config["RSeQC-check_strandedness_assessable_from_infer_experiment"])


        with vp.component_start(
            name="RSEM Counts",
            description="Dataset wide checks including outliers detection",
            ):
            with vp.payload(payloads=[
                {"dataset":dataset, "sample_component": "geneCounts"}
            ]) as ADD:
                ADD(check_for_outliers, config=config["RSEM Counts-check_for_outliers"])

        sample: BulkRNASeqSample
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
                            config=config[
                                "Raw Reads By Sample-check_fastqgz_file_contents"
                            ],
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
                            config=config[
                                "Trim Reads By Sample-check_fastqgz_file_contents"
                            ],
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
                            },
                        ]
                    ) as ADD:
                        ADD(
                            check_thresholds,
                            config=config[
                                "STAR Alignments By Sample-check_thresholds-Mapped"
                            ],
                            description="Check that mapping rates are reasonable, specifically most reads map to the target genome",
                        )
                        ADD(
                            check_thresholds,
                            config=config[
                                "STAR Alignments By Sample-check_thresholds-MultiMapped"
                            ],
                            description="Check that mapping rates are reasonable, specifically that a considerable amount of reads multimap to the target genome",
                        )

                with vp.component_start(
                    name="RSeQC By Sample",
                    description="RNASeq QA outputs",
                ):
                    with vp.component_start(
                        name="geneBody_coverage",
                        description="Assess integrity of transcripts and library prep signatures",
                    ):
                        with vp.payload(
                            payloads=[
                                {"file": lambda: sample.rSeQCAnalysis.geneBodyCoverageMultiQCDirZIP.path},
                            ]
                        ) as ADD:
                            ADD(check_file_exists)
                        with vp.payload(
                            payloads=[
                                {"input_dir": lambda: sample.rSeQCAnalysis.geneBodyCoverageOut.path},
                            ]
                        ) as ADD:
                            ADD(check_genebody_coverage_output)

                    with vp.component_start(
                        name="infer_experiment",
                        description="Assess strandedness of transcripts based on gene annotations",
                    ):
                        with vp.payload(
                            payloads=[
                                {"file": lambda: sample.rSeQCAnalysis.inferExperimentMultiQCDirZIP.path},
                                {"file": lambda: sample.rSeQCAnalysis.inferExperimentOut.path},
                            ]
                        ) as ADD:
                            ADD(check_file_exists)

                    with vp.component_start(
                        name="inner_distance",
                        description="Reports on distance between mate reads based on gene annotations",
                        skip=(not dataset.metadata.paired_end)
                    ):
                        with vp.payload(
                            payloads=[
                                {"file": lambda: sample.rSeQCAnalysis.innerDistanceMultiQCDirZIP.path},
                            ]
                        ) as ADD:
                            ADD(check_file_exists)
                        with vp.payload(
                            payloads=[
                                {"input_dir": lambda: sample.rSeQCAnalysis.innerDistanceOut.path},
                            ]
                        ) as ADD:
                            ADD(check_inner_distance_output)


                    with vp.component_start(
                        name="read_distribution",
                        description="Assess average element makeup of transcript",
                    ):
                        with vp.payload(
                            payloads=[
                                {"file": lambda: sample.rSeQCAnalysis.readDistributionMultiQCDirZIP.path},
                                {"file": lambda: sample.rSeQCAnalysis.readDistributionOut.path},
                            ]
                        ) as ADD:
                            ADD(check_file_exists)
    # return protocol object without running or generating a report
    if defer_run:
        return vp

    vp.run()

    # return report
    return vp.report(**report_args)

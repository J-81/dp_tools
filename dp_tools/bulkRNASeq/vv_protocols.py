from collections import defaultdict
import enum
import os
from pathlib import Path
from typing import Dict, List, Set, Tuple, Union
import yaml
import logging

from dp_tools.core.entity_model2 import Dataset

log = logging.getLogger(__name__)

from dp_tools.core.check_model import ValidationProtocol
from dp_tools.bulkRNASeq.checks import *  # normally this isn't ideal, however, as a library of functions this seems reasonable

CONFIG = {
    "Metadata-check_metadata_attributes_exist": {
        "expected_attrs": ["paired_end", "has_ercc"]
    },
    "Raw Reads-check_for_outliers": {
        "mqc_module": "FastQC",
        "mqc_plot": "general_stats",
        "mqc_keys": [
            "percent_gc",
            "avg_sequence_length",
            "total_sequences",
            "percent_duplicates",
        ],
        "thresholds": [
            {"code": "YELLOW", "stdev_threshold": 2, "middle_fcn": "median"},
            {"code": "RED", "stdev_threshold": 4, "middle_fcn": "median"},
        ],
    },
    "Trim Reads-check_for_outliers": {
        "mqc_module": "FastQC",
        "mqc_plot": "general_stats",
        "mqc_keys": [
            "percent_gc",
            "avg_sequence_length",
            "total_sequences",
            "percent_duplicates",
        ],
        "thresholds": [
            {"code": "YELLOW", "stdev_threshold": 2, "middle_fcn": "median"},
            {"code": "RED", "stdev_threshold": 4, "middle_fcn": "median"},
        ],
    },
    "Raw Reads By Sample-check_fastqgz_file_contents": {
        "count_lines_to_check": 200000000
    },
    "Trim Reads By Sample-check_fastqgz_file_contents": {
        "count_lines_to_check": 200000000
    },
    "STAR Alignments By Sample-check_thresholds-Mapped": {
        "mqc_key": "STAR",
        "stat_string": "uniquely_mapped_percent + multimapped_percent",
        "thresholds": [
            {"code": "YELLOW", "type": "lower", "value": 70},
            {"code": "RED", "type": "lower", "value": 50},
        ],
    },
    "STAR Alignments By Sample-check_thresholds-MultiMapped": {
        "mqc_key": "STAR",
        "stat_string": "multimapped_toomany_percent + multimapped_percent",
        "thresholds": [
            {"code": "YELLOW", "type": "lower", "value": 30},
            {"code": "RED", "type": "lower", "value": 15},
        ],
    },
    "STAR Alignments-check_for_outliers": {
        "mqc_module": "STAR",
        "mqc_plot": "general_stats",
        "mqc_keys": [
            "uniquely_mapped_percent",
            "avg_mapped_read_length",
            "mismatch_rate",
            "deletion_rate",
            "deletion_length",
            "insertion_rate",
            "insertion_length",
            "multimapped_percent",
            "multimapped_toomany_percent",
            "unmapped_mismatches_percent",
            "unmapped_tooshort_percent",
            "unmapped_other_percent",
        ],
        "thresholds": [
            {"code": "YELLOW", "stdev_threshold": 2, "middle_fcn": "median"},
            {"code": "RED", "stdev_threshold": 4, "middle_fcn": "median"},
        ],
    },
    "RSeQC-check_for_outliers-geneBody_coverage": {
        "mqc_module": "RSeQC",
        "mqc_plot": "Gene Body Coverage",
        "mqc_keys": ["_ALL"],
        "thresholds": [
            {"code": "YELLOW", "stdev_threshold": 2, "middle_fcn": "median"},
            {"code": "RED", "stdev_threshold": 4, "middle_fcn": "median"},
        ],
    },
    "RSeQC-check_for_outliers-infer_experiment": {
        "mqc_module": "RSeQC",
        "mqc_plot": "Infer experiment",
        "mqc_keys": ["_ALL"],
        "thresholds": [
            {"code": "YELLOW", "stdev_threshold": 2, "middle_fcn": "median"},
            {"code": "RED", "stdev_threshold": 4, "middle_fcn": "median"},
        ],
    },
    "RSeQC-check_for_outliers-inner_distance": {
        "mqc_module": "RSeQC",
        "mqc_plot": "Inner Distance",
        "mqc_keys": ["_ALL"],
        "thresholds": [
            {"code": "YELLOW", "stdev_threshold": 2, "middle_fcn": "median"},
            {"code": "RED", "stdev_threshold": 4, "middle_fcn": "median"},
        ],
    },
    "RSeQC-check_for_outliers-read_distribution": {
        "mqc_module": "RSeQC",
        "mqc_plot": "Read Distribution",
        "mqc_keys": ["_ALL"],
        "thresholds": [
            {"code": "YELLOW", "stdev_threshold": 2, "middle_fcn": "median"},
            {"code": "RED", "stdev_threshold": 4, "middle_fcn": "median"},
        ],
    },
    "RSeQC-check_strandedness_assessable_from_infer_experiment": {
        "stranded_assessment_range": {"max": 100, "min": 75},
        "unstranded_assessment_range": {"min": 40, "max": 60},
        "valid_dominant_strandedness_assessments": [
            "Sense (% Tags)",
            "Antisense (% Tags)",
        ],
    },
    "RSEM Counts-check_for_outliers": {
        "mqc_module": "Rsem",
        "mqc_plot": "general_stats",
        "mqc_keys": [
            "Unalignable",
            "Alignable",
            "Filtered",
            "Total",
            "alignable_percent",
            "Unique",
            "Multi",
            "Uncertain",
        ],
        "thresholds": [
            {"code": "YELLOW", "stdev_threshold": 2, "middle_fcn": "median"},
            {"code": "RED", "stdev_threshold": 4, "middle_fcn": "median"},
        ],
    },
}


def validate_bulkRNASeq(
    dataset: Dataset,
    config_path: Path = None,
    run_args: dict = None,
    report_args: dict = None,
    protocol_args: dict = None,
    defer_run: bool = False,
) -> Union[ValidationProtocol, ValidationProtocol.Report]:

    if config_path is not None:
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)
    else:
        config = CONFIG

    if run_args is None:
        run_args = dict()

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
            with vp.payload(payloads=[{"dataset": dataset}]):
                vp.add(
                    check_metadata_attributes_exist,
                    config=config["Metadata-check_metadata_attributes_exist"],
                )

        with vp.component_start(
            name="Raw Reads", description="Raw Reads Outliers Detection"
        ):
            with vp.payload(
                payloads=[{"dataset": dataset, "sample_component": "rawReads"}]
                if not dataset.metadata['paired_end']
                else [
                    {"dataset": dataset, "sample_component": "rawForwardReads"},
                    {"dataset": dataset, "sample_component": "rawReverseReads"},
                ]
            ):
                vp.add(
                    check_for_outliers, config=config["Raw Reads-check_for_outliers"]
                )

        with vp.component_start(
            name="Trim Reads", description="Trimmed Reads Outliers Detection"
        ):
            with vp.payload(
                payloads=[{"dataset": dataset, "sample_component": "trimReads"}]
                if not dataset.metadata['paired_end']
                else [
                    {"dataset": dataset, "sample_component": "trimForwardReads"},
                    {"dataset": dataset, "sample_component": "trimReverseReads"},
                ]
            ):
                vp.add(
                    check_for_outliers,
                    config=config["Trim Reads-check_for_outliers"],
                )

        with vp.component_start(
            name="STAR Alignments",
            description="Dataset wide checks including outliers detection",
        ):
            with vp.payload(
                payloads=[{"dataset": dataset, "sample_component": "genomeAlignments"}]
            ):
                vp.add(
                    check_for_outliers,
                    config=config["STAR Alignments-check_for_outliers"],
                )

        with vp.component_start(
            name="RSeQC",
            description="RSeQC submodule outliers checking and other submodule specific dataset wide checks",
        ):
            with vp.payload(
                payloads=[{"dataset": dataset, "sample_component": "rSeQCAnalysis"}]
            ):
                vp.add(
                    check_for_outliers,
                    description="Check for outliers in geneBody Coverage",
                    config=config["RSeQC-check_for_outliers-geneBody_coverage"],
                )
                vp.add(
                    check_for_outliers,
                    description="Check for outliers in infer experiment",
                    config=config["RSeQC-check_for_outliers-infer_experiment"],
                )
                vp.add(
                    check_for_outliers,
                    description="Check for outliers in inner distance",
                    config=config["RSeQC-check_for_outliers-inner_distance"],
                    skip=(not dataset.metadata['paired_end']),
                )
                vp.add(
                    check_for_outliers,
                    description="Check for outliers in read distribution",
                    config=config["RSeQC-check_for_outliers-read_distribution"],
                )

            with vp.payload(payloads=[{"dataset": dataset}]):
                vp.add(
                    check_strandedness_assessable_from_infer_experiment,
                    config=config[
                        "RSeQC-check_strandedness_assessable_from_infer_experiment"
                    ],
                )

        with vp.component_start(
            name="RSEM Counts",
            description="Dataset wide checks including outliers detection",
        ):
            with vp.payload(
                payloads=[{"dataset": dataset, "sample_component": "geneCounts"}]
            ):
                vp.add(
                    check_for_outliers, config=config["RSEM Counts-check_for_outliers"]
                )

        with vp.component_start(
            name="Unnormalized Gene Counts",
            description="Validate normalization related output",
        ):

            with vp.payload(
                payloads=[
                    {
                        "unnormalizedCountTable": lambda: dataset.starGeneCounts.unnormalizedCounts.path,
                        "samplewise_tables": lambda: {
                            s.name: s.genomeAlignments.readsPerGeneTable.path
                            for s in dataset.samples.values()
                        },
                    },
                ]
            ):
                vp.add(
                    check_aggregate_star_unnormalized_counts_table_values_against_samplewise_tables
                )
            with vp.payload(
                payloads=[
                    {
                        "unnormalizedCountTable": lambda: dataset.rsemGeneCounts.unnormalizedCounts.path,
                        "samplewise_tables": lambda: {
                            s.name: s.geneCounts.genesResults.path
                            for s in dataset.samples.values()
                        },
                    },
                ]
            ):
                vp.add(
                    check_aggregate_rsem_unnormalized_counts_table_values_against_samplewise_tables
                )
                vp.add(check_ERCC_subgroup_representation, skip=(not dataset.metadata['has_ERCC']))

        with vp.component_start(
            name="DGE Metadata",
            description="",
            ):
            
            with vp.component_start(
                name="Sample Table",
                description="",
                ):
                with vp.payload(payloads=[
                    {
                    'runsheet': lambda: dataset.data_assets['runsheet'],
                    'sampleTable': lambda: dataset.data_assets['sample table']
                    }
                ]):
                    vp.add(check_sample_table_against_runsheet, config={"all_samples_required": True})
                    vp.add(check_sample_table_for_correct_group_assignments)
                    
            with vp.component_start(
                name="Contrasts Tables",
                description="",
                ):
                with vp.payload(payloads=[
                    {
                    'runsheet': lambda: dataset.data_assets['runsheet'],
                    'contrasts_table': lambda: dataset.data_assets['DESeq2 contrasts table']
                    }
                ]):
                    vp.add(check_contrasts_table_headers)
                    vp.add(check_contrasts_table_rows)

        with vp.component_start(
            name="DGE Metadata ERCC",
            description="",
            skip=(not dataset.metadata['has_ERCC'])
            ):

            with vp.component_start(
                name="Sample Table",
                description="",
                ):
                with vp.payload(payloads=[
                    {
                    'runsheet': lambda: dataset.data_assets['runsheet'],
                    'sampleTable': lambda: dataset.data_assets['ERCC sample table']
                    }
                ]):
                    vp.add(check_sample_table_against_runsheet, config={"all_samples_required": False})
                    vp.add(check_sample_table_for_correct_group_assignments) 

            with vp.component_start(
                name="Contrasts Tables",
                description="",
                ):
                with vp.payload(payloads=[
                    {
                    'runsheet': lambda: dataset.data_assets['runsheet'],
                    'contrasts_table': lambda: dataset.data_assets['ERCC normalized DESeq2 contrasts table']
                    }
                ]):
                    vp.add(check_contrasts_table_headers)
                    vp.add(check_contrasts_table_rows)
                    
        with vp.component_start(
            name="DGE Output",
            description="",
            ):
            with vp.payload(
                payloads=[
                    {
                        "rsem_table_path": lambda: dataset.data_assets['rsem unnormalized counts table'],
                        "deseq2_table_path": lambda: dataset.data_assets['DESeq2 unnormalized counts table'],
                    }
                ]
            ):
                vp.add(check_rsem_counts_and_unnormalized_tables_parity)

            with vp.payload(payloads=[
                {
                'organism': lambda: dataset.metadata['organism'],
                'samples': lambda: set(dataset.samples),
                'dge_table': lambda: dataset.data_assets['DESeq2 annotated DGE table'],
                'runsheet': lambda: dataset.data_assets['runsheet'],
                }
            ]):
                vp.add(check_dge_table_annotation_columns_exist)
                vp.add(check_dge_table_sample_columns_exist)
                vp.add(check_dge_table_sample_columns_constraints)
                vp.add(check_dge_table_group_columns_exist)
                vp.add(check_dge_table_group_columns_constraints)
                vp.add(check_dge_table_comparison_statistical_columns_exist)
                vp.add(check_dge_table_group_statistical_columns_constraints)
                vp.add(check_dge_table_fixed_statistical_columns_exist)
                vp.add(check_dge_table_fixed_statistical_columns_constraints)
                vp.add(check_dge_table_log2fc_within_reason)

            with vp.component_start(
                name="Viz Tables",
                description="Extended from the dge tables",
                ):
                with vp.payload(payloads=[
                    {
                    'organism': lambda: dataset.metadata['organism'],
                    'samples': lambda: set(dataset.samples),
                    'dge_table': lambda: dataset.data_assets['DESeq2 annotated DGE extended for viz table'],
                    'runsheet': lambda: dataset.data_assets['runsheet'],
                    }
                ]):
                    vp.add(check_dge_table_annotation_columns_exist)
                    vp.add(check_dge_table_sample_columns_exist)
                    vp.add(check_dge_table_sample_columns_constraints)
                    vp.add(check_dge_table_group_columns_exist)
                    vp.add(check_dge_table_group_columns_constraints)
                    vp.add(check_dge_table_comparison_statistical_columns_exist)
                    vp.add(check_dge_table_group_statistical_columns_constraints)
                    vp.add(check_dge_table_fixed_statistical_columns_exist)
                    vp.add(check_dge_table_fixed_statistical_columns_constraints)
                    vp.add(check_dge_table_log2fc_within_reason)
                    vp.add(check_viz_table_columns_exist)
                    vp.add(check_viz_table_columns_constraints)

                with vp.payload(payloads=[
                    {
                    'samples': lambda: set(dataset.samples),
                    'pca_table': lambda: dataset.data_assets['DESeq2 viz PCA table'],
                    }
                ]):
                    vp.add(check_viz_pca_table_index_and_columns_exist)

        with vp.component_start(
            name="DGE Output ERCC",
            description="",
            skip=(not dataset.metadata['has_ERCC'])
            ):
            with vp.payload(payloads=[
                {
                'organism': lambda: dataset.metadata['organism'],
                'samples': lambda: set(pd.read_csv(dataset.data_assets['ERCC sample table'], index_col=0).index),
                'dge_table': lambda: dataset.data_assets['ERCC normalized DESeq2 annotated DGE table'],
                'runsheet': lambda: dataset.data_assets['runsheet'],
                }
            ]):
                vp.add(check_dge_table_annotation_columns_exist)
                vp.add(check_dge_table_sample_columns_exist)
                vp.add(check_dge_table_sample_columns_constraints)
                vp.add(check_dge_table_group_columns_exist)
                vp.add(check_dge_table_group_columns_constraints)
                vp.add(check_dge_table_comparison_statistical_columns_exist)
                vp.add(check_dge_table_group_statistical_columns_constraints)
                vp.add(check_dge_table_fixed_statistical_columns_exist)
                vp.add(check_dge_table_fixed_statistical_columns_constraints)
                vp.add(check_dge_table_log2fc_within_reason)

            with vp.component_start(
                name="Viz Tables",
                description="Extended from the dge tables",
                ):
                with vp.payload(payloads=[
                    {
                    'organism': lambda: dataset.metadata['organism'],
                    'samples': lambda: set(pd.read_csv(dataset.normalizedGeneCounts.erccSampleTableCSV.path, index_col=0).index),
                    'dge_table': lambda: dataset.differentialGeneExpressionERCC.visualizationTableCSV.path,
                    'runsheet': lambda: dataset.metadata.runsheet.path
                    }
                ]):
                    vp.add(check_dge_table_annotation_columns_exist)
                    vp.add(check_dge_table_sample_columns_exist)
                    vp.add(check_dge_table_sample_columns_constraints)
                    vp.add(check_dge_table_group_columns_exist)
                    vp.add(check_dge_table_group_columns_constraints)
                    vp.add(check_dge_table_comparison_statistical_columns_exist)
                    vp.add(check_dge_table_group_statistical_columns_constraints)
                    vp.add(check_dge_table_fixed_statistical_columns_exist)
                    vp.add(check_dge_table_fixed_statistical_columns_constraints)
                    vp.add(check_dge_table_log2fc_within_reason)
                    vp.add(check_viz_table_columns_exist)
                    vp.add(check_viz_table_columns_constraints)

                with vp.payload(payloads=[
                    {
                    'samples': lambda: set(pd.read_csv(dataset.normalizedGeneCounts.erccSampleTableCSV.path, index_col=0).index),
                    'pca_table': lambda: dataset.differentialGeneExpressionERCC.visualizationPCATableCSV.path,
                    }
                ]):
                    vp.add(check_viz_pca_table_index_and_columns_exist)

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
                            if dataset.metadata['paired_end']
                            else [{"file": lambda: sample.rawReads.fastqGZ.path}]
                        )
                    ):
                        vp.add(check_file_exists, description="Check reads files exist")
                        vp.add(
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
                    ):
                        vp.add(
                            check_forward_and_reverse_reads_counts_match,
                            skip=(not dataset.metadata['paired_end']),
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
                                if dataset.metadata['paired_end']
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
                        ):
                            vp.add(check_file_exists)

                with vp.component_start(
                    name="Trimmed Reads By Sample", description="Trimmed reads"
                ):
                    with vp.payload(
                        payloads=(
                            [
                                {"file": lambda: sample.data_assets.get("trimmed forward reads fastq GZ")},
                                {"file": lambda:  sample.data_assets.get("trimmed reverse reads fastq GZ")},
                            ]
                            if dataset.metadata['paired_end']
                            else [{"file": lambda:  sample.data_assets.get("trimmed reads fastq GZ")}]
                        )
                    ):
                        vp.add(check_file_exists, description="Check reads files exist")
                        vp.add(
                            check_fastqgz_file_contents,
                            config=config[
                                "Trim Reads By Sample-check_fastqgz_file_contents"
                            ],
                        )

                    with vp.payload(
                        payloads=[
                            {
                                "fwd_reads": lambda: sample.trimForwardReads,
                                "rev_reads": lambda: sample.trimReverseReads,
                            },
                        ],
                    ):
                        # vp.add(
                        #     check_forward_and_reverse_reads_counts_match,
                        #     skip=(not dataset.metadata['paired_end']),
                        # )
                        ... # TODO: reimplement

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
                                if dataset.metadata['paired_end']
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
                        ):
                            vp.add(check_file_exists)

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
                                if dataset.metadata['paired_end']
                                else [
                                    {
                                        "file": lambda: sample.trimReads.trimmingReportTXT.path
                                    },
                                ]
                            )
                        ):
                            vp.add(
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
                    ):
                        vp.add(
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
                    ):
                        # vp.add(check_bam_file_integrity, config={"samtools_bin":#TODO: fill with private config file})
                        vp.add(
                            check_bam_file_integrity,
                            config={"samtools_bin": 'samtools'}, # assumes accessible on path already
                        )

                    with vp.payload(
                        payloads=[
                            {
                                "component": lambda: sample.genomeAlignments,
                            },
                        ]
                    ):
                        vp.add(
                            check_thresholds,
                            config=config[
                                "STAR Alignments By Sample-check_thresholds-Mapped"
                            ],
                            description="Check that mapping rates are reasonable, specifically most reads map to the target genome",
                        )
                        vp.add(
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
                                {
                                    "file": lambda: sample.rSeQCAnalysis.geneBodyCoverageMultiQCDirZIP.path
                                },
                            ]
                        ):
                            vp.add(check_file_exists)
                        with vp.payload(
                            payloads=[
                                {
                                    "input_dir": lambda: sample.rSeQCAnalysis.geneBodyCoverageOut.path
                                },
                            ]
                        ):
                            vp.add(check_genebody_coverage_output)

                    with vp.component_start(
                        name="infer_experiment",
                        description="Assess strandedness of transcripts based on gene annotations",
                    ):
                        with vp.payload(
                            payloads=[
                                {
                                    "file": lambda: sample.rSeQCAnalysis.inferExperimentMultiQCDirZIP.path
                                },
                                {
                                    "file": lambda: sample.rSeQCAnalysis.inferExperimentOut.path
                                },
                            ]
                        ):
                            vp.add(check_file_exists)

                    with vp.component_start(
                        name="inner_distance",
                        description="Reports on distance between mate reads based on gene annotations",
                        skip=(not dataset.metadata['paired_end']),
                    ):
                        with vp.payload(
                            payloads=[
                                {
                                    "file": lambda: sample.rSeQCAnalysis.innerDistanceMultiQCDirZIP.path
                                },
                            ]
                        ):
                            vp.add(check_file_exists)
                        with vp.payload(
                            payloads=[
                                {
                                    "input_dir": lambda: sample.rSeQCAnalysis.innerDistanceOut.path
                                },
                            ]
                        ):
                            vp.add(check_inner_distance_output)

                    with vp.component_start(
                        name="read_distribution",
                        description="Assess average element makeup of transcript",
                    ):
                        with vp.payload(
                            payloads=[
                                {
                                    "file": lambda: sample.rSeQCAnalysis.readDistributionMultiQCDirZIP.path
                                },
                                {
                                    "file": lambda: sample.rSeQCAnalysis.readDistributionOut.path
                                },
                            ]
                        ):
                            vp.add(check_file_exists)
    # return protocol object without running or generating a report
    if defer_run:
        return vp

    vp.run(**run_args)

    # return report
    return vp.report(**report_args)

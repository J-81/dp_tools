import os
from pathlib import Path
from dp_tools.bulkRNASeq.checks import (
    check_ERCC_subgroup_representation,
    check_bam_file_integrity,
    check_forward_and_reverse_reads_counts_match,
    check_rsem_counts_and_unnormalized_tables_parity,
)
from dp_tools.bulkRNASeq.entity import BulkRNASeqDataset
from dp_tools.core.check_model import FlagCode


def test_check_rsem_counts_and_unnormalized_tables_parity(glds48_dataSystem_STAGE04):
    dataset = glds48_dataSystem_STAGE04.dataset
    res = check_rsem_counts_and_unnormalized_tables_parity(
        rsem_table_path=dataset.rsemGeneCounts.unnormalizedCounts.path,
        deseq2_table_path=dataset.normalizedGeneCounts.unnormalizedCountsCSV.path,
    )

    assert res["code"]
    assert res["message"]


def test_check_forward_and_reverse_reads_counts_match(glds194_dataSystem_STAGE00):
    dataset = glds194_dataSystem_STAGE00.dataset

    for sample in dataset.samples.values():
        res = check_forward_and_reverse_reads_counts_match(
            fwd_reads=sample.rawForwardReads, rev_reads=sample.rawReverseReads
        )
        assert res["code"]
        assert res["message"]


def test_check_bam_file_integrity(glds48_dataSystem_STAGE04):
    dataset = glds48_dataSystem_STAGE04.dataset

    for sample in dataset.samples.values():
        res = check_bam_file_integrity(
            file=sample.genomeAlignments.alignedToTranscriptomeBam.path,
            samtools_bin=Path(os.environ.get("SAMTOOLS_BIN")),
        )
        assert res["code"] == FlagCode.GREEN
        assert res["message"]

        res = check_bam_file_integrity(
            file=sample.genomeAlignments.logFinal.path,
            samtools_bin=Path(os.environ.get("SAMTOOLS_BIN")),
        )

        assert res["code"] == FlagCode.HALT
        assert res["message"]


def test_check_ERCC_group_represention(glds194_dataSystem_STAGE04):
    dataset: BulkRNASeqDataset = glds194_dataSystem_STAGE04.dataset

    res = check_ERCC_subgroup_representation(
        unnormalizedCountTable=dataset.rsemGeneCounts.unnormalizedCounts.path
    )

    assert res["code"] == FlagCode.RED
    assert res["message"]

import os
from pathlib import Path
from dp_tools.bulkRNASeq.checks import (
    check_bam_file_integrity,
    check_forward_and_reverse_reads_counts_match,
)
from dp_tools.core.check_model import FlagCode


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

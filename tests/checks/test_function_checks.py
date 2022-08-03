import os
from pathlib import Path
from dp_tools.bulkRNASeq.checks import *
from dp_tools.core.check_model import FlagCode


def test_check_rsem_counts_and_unnormalized_tables_parity(glds48_dataSystem):
    dataset = glds48_dataSystem.dataset
    res = check_rsem_counts_and_unnormalized_tables_parity(
        rsem_table_path=dataset.data_assets["rsem unnormalized counts table"].path,
        deseq2_table_path=dataset.data_assets["DESeq2 unnormalized counts table"].path,
    )

    assert (
        res["code"] == FlagCode.HALT
    )  # The test data was truncate processed and counts won't match
    assert res["message"]


def test_check_forward_and_reverse_reads_counts_match(glds194_dataSystem):
    dataset = glds194_dataSystem.dataset

    for sample in dataset.samples.values():
        res = check_forward_and_reverse_reads_counts_match(
            sample=sample,
            reads_key_1="raw forward reads fastQC ZIP",
            reads_key_2="raw reverse reads fastQC ZIP",
        )
        assert res["code"] == FlagCode.GREEN
        assert res["message"]


def test_check_bam_file_integrity(glds48_dataSystem):
    dataset = glds48_dataSystem.dataset

    for sample in dataset.samples.values():
        res = check_bam_file_integrity(
            file=sample.data_assets["aligned SortedByCoord Bam"].path,
        )
        assert res["code"] == FlagCode.GREEN
        assert res["message"]

        res = check_bam_file_integrity(
            file=sample.data_assets["aligned log Full"].path,
        )

        assert res["code"] == FlagCode.HALT
        assert res["message"]


def test_check_ERCC_group_represention(glds194_dataSystem):
    dataset = glds194_dataSystem.dataset

    res = check_ERCC_subgroup_representation(
        unnormalizedCountTable=dataset.data_assets[
            "rsem unnormalized counts table"
        ].path
    )

    assert res["code"] == FlagCode.RED
    assert res["message"]  # TODO: Fix message


def test_check_dge_table_log2fc_within_reason(glds194_dataSystem):
    dataset = glds194_dataSystem.dataset

    res = check_dge_table_log2fc_within_reason(
        dge_table=dataset.data_assets["DESeq2 annotated DGE table"].path,
        runsheet=dataset.data_assets["runsheet"].path,
    )

    assert res["code"] == FlagCode.GREEN
    assert res["message"]


def test_check_sample_table_against_runsheet(glds194_dataSystem):
    dataset = glds194_dataSystem.dataset

    res = check_sample_table_against_runsheet(
        runsheet=dataset.data_assets["runsheet"].path,
        sampleTable=dataset.data_assets["sample table"].path,
        all_samples_required=True,
    )

    assert res["code"] == FlagCode.GREEN
    assert res["message"]

    res = check_sample_table_against_runsheet(
        runsheet=dataset.data_assets["runsheet"].path,
        sampleTable=dataset.data_assets["ERCC sample table"].path,
        all_samples_required=False,
    )

    assert res["code"] == FlagCode.GREEN
    assert res["message"]


def test_check_contrasts_table_rows(glds48_dataSystem):
    dataset = glds48_dataSystem.dataset

    res = check_contrasts_table_rows(
        contrasts_table=dataset.data_assets["DESeq2 contrasts table"].path
    )

    assert res["code"] == FlagCode.GREEN
    assert res["message"]


def test_check_sample_table_for_correct_group_assignments(glds194_dataSystem):
    dataset = glds194_dataSystem.dataset

    res = check_sample_table_for_correct_group_assignments(
        runsheet=dataset.data_assets["runsheet"].path,
        sampleTable=dataset.data_assets["sample table"].path,
    )

    assert res["code"] == FlagCode.GREEN
    assert res["message"]

    res = check_sample_table_for_correct_group_assignments(
        runsheet=dataset.data_assets["runsheet"].path,
        sampleTable=dataset.data_assets["ERCC sample table"].path,
    )

    assert res["code"] == FlagCode.GREEN
    assert res["message"]


def test_check_dge_table_group_columns_constraints(glds194_dataSystem):
    dataset = glds194_dataSystem.dataset

    payload = {
        "organism": dataset.metadata["organism"],
        "samples": set(dataset.samples),
        "dge_table": dataset.data_assets["DESeq2 annotated DGE table"].path,
        "runsheet": dataset.data_assets["runsheet"].path,
    }

    res = check_dge_table_group_columns_constraints(**payload)

    assert res["code"] == FlagCode.GREEN
    assert res["message"]

    payload = {
        "organism": dataset.metadata["organism"],
        "samples": set(
            pd.read_csv(
                dataset.data_assets["ERCC sample table"].path, index_col=0
            ).index
        ),
        "dge_table": dataset.data_assets["DESeq2 annotated DGE table"].path,
        "runsheet": dataset.data_assets['runsheet'].path,
    }

    res = check_dge_table_group_columns_constraints(**payload)

    assert res["code"] == FlagCode.GREEN
    assert res["message"]

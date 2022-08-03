""" Tests for validation report results, relies on test for loaders passing """
from decimal import DivisionByZero
from pathlib import Path
import os

import pytest
from dp_tools.bulkRNASeq.vv_protocols import (
    validate_bulkRNASeq,
)
from dp_tools.core.check_model import ValidationProtocol
from dp_tools.core.loaders import load_data


@pytest.fixture(autouse=True)
def mock_dev_exceptions(monkeypatch):
    monkeypatch.setattr(
        "dp_tools.core.check_model.ALLOWED_DEV_EXCEPTIONS", (DivisionByZero)
    )  # ensure unhandled developer exceptions are raised


def pseudo_fingerprint(df):
    return (
        df["code"].apply(lambda v: v.value).sum()
        + df["code"].apply(lambda v: v.value).mean()
        + df.shape[0] * df.shape[1]
    )


# TODO: Part of an alternative framework not fully implemented

import logging

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)


@pytest.mark.parametrize(
    "data_asset_keys,components,expected_flag_table_shape,expected_outlier_table_shape,expected_flag_table_fingerprint,expected_halt_flag_count",
    [
        pytest.param(
            ["demuliplexed paired end raw data", "qc reports for paired end raw data"],
            ["Metadata", "Raw Reads By Sample", "Raw Reads"],
            (124, 7),
            (4, 1),
            3418.4032258064517,
            0,
            id="Raw Reads Checks Only",
        ),
        pytest.param(
            [
                "paired end trimmed reads",
                "qc reports for paired end trimmed reads data",
            ],
            ["Trimmed Reads By Sample", "Trim Reads"],
            (173, 7),
            (9, 1),
            4741.289017341041,
            0,
            id="Trimmed Reads Checks Only",
        ),
        pytest.param(
            ["STAR alignments"],
            ["STAR Alignments", "STAR Alignments By Sample"],
            (175, 7),
            (7, 1),
            5569.685714285714,
            0,
            id="STAR Alignments Checks Only",
        ),
        pytest.param(
            ["RSeQC output for paired end data"],
            ["RSeQC", "RSeQC By Sample"],
            (92, 7),
            (13, 1),
            2635.413043478261,
            1,  # RSeQC strandedness is ambiguous
            id="RSeQC Checks Only",
        ),
        pytest.param(
            ["RSEM counts"],
            ["RSEM Counts"],
            (45, 7),
            (3, 1),
            1275.888888888889,
            0,
            id="RSEM Counts Checks Only",
        ),
        pytest.param(
            ["RSEM counts", "STAR alignments"],
            ["Unnormalized Gene Counts"],
            (168, 7),
            (0, 0),
            4616.357142857143,
            0,
            id="Unnormalized Gene Counts Checks Only",
        ),
        pytest.param(
            ["DGE Output", "ERCC DGE Output", "RSEM Output"],
            ["DGE Metadata", "DGE Metadata ERCC", "DGE Output", "DGE Output ERCC"],
            (73, 7),
            (0, 0),
            2082.232876712329,
            1,  # RSEM and DGE parity won't match due to dummy counts in DGE output test data
            id="DGE Checks Only",
        ),
        pytest.param(
            ["is paired end full", "ERCC DGE Output"],
            None,  # This evaluates to meaning running all components
            (679, 7),
            (34, 5),
            19414.561119293077,
            2,  # RSEM and DGE parity won't match due to dummy counts in DGE output test data / RSeQC strandedness is ambiguous
            id="Run all checks",
        ),
    ],
)
def test_updated_protocol_model_paired_end(
    glds194_test_dir,
    data_asset_keys,
    components,
    expected_flag_table_shape,
    expected_outlier_table_shape,
    expected_flag_table_fingerprint,
    expected_halt_flag_count,
):
    datasystem = load_data(
        key_sets=data_asset_keys,
        config=("bulkRNASeq", "Latest"),
        root_path=(glds194_test_dir),
        runsheet_path=(
            glds194_test_dir / "Metadata/GLDS-194_bulkRNASeq_v1_runsheet.csv"
        ),
    )

    report = validate_bulkRNASeq(
        datasystem.dataset,
        report_args={"include_skipped": False},
        protocol_args={"run_components": components},
    )

    assert (
        sum(report["flag_table"]["code_level"] >= 80) == expected_halt_flag_count
    ), f"Found more than expected HALT+ flags: {report['flag_table'].loc[report['flag_table']['code_level'] >= 80].to_dict()}"

    assert (
        report["flag_table"].shape,
        report["outliers"].shape,
        pseudo_fingerprint(report["flag_table"]),
    ) == (
        expected_flag_table_shape,
        expected_outlier_table_shape,
        expected_flag_table_fingerprint,
    )


@pytest.mark.parametrize(
    "data_asset_keys,components,expected_flag_table_shape,expected_outlier_table_shape,expected_flag_table_fingerprint,expected_halt_flag_count",
    [
        pytest.param(
            ["demuliplexed single end raw data", "qc reports for single end raw data"],
            ["Metadata", "Raw Reads By Sample", "Raw Reads"],
            (61, 7),
            (1, 1),
            1677.1639344262296,
            0,
            id="Raw Reads Checks Only",
        ),
        pytest.param(
            [
                "single end trimmed reads",
                "qc reports for single end trimmed reads data",
            ],
            ["Trimmed Reads By Sample", "Trim Reads"],
            (86, 7),
            (2, 1),
            2352.116279069767,
            0,
            id="Trimmed Reads Checks Only",
        ),
        pytest.param(
            ["STAR alignments"],
            ["STAR Alignments", "STAR Alignments By Sample"],
            (187, 7),
            (7, 1),
            5923.545454545455,
            0,
            id="STAR Alignments Checks Only",
        ),
        pytest.param(
            ["RSeQC output for single end data"],
            ["RSeQC", "RSeQC By Sample"],
            (66, 7),
            (12, 1),
            1862.909090909091,
            0,
            id="RSeQC Checks Only",
        ),
        pytest.param(
            ["RSEM counts"],
            ["RSEM Counts"],
            (47, 7),
            (1, 1),
            1299.2127659574467,
            0,
            id="RSEM Counts Checks Only",
        ),
        pytest.param(
            ["RSEM counts", "STAR alignments"],
            ["Unnormalized Gene Counts"],
            (178, 7),
            (0, 0),
            4886.337078651685,
            0,
            id="Unnormalized Gene Counts Checks Only",
        ),
        pytest.param(
            ["DGE Output", "RSEM Output"],
            ["DGE Metadata", "DGE Metadata ERCC", "DGE Output", "DGE Output ERCC"],
            (39, 7),
            (0, 0),
            1134.5384615384614,
            1,  # test data does not have parity between DGE counts and RSEM counts
            id="DGE Checks Only",
        ),
        pytest.param(
            ["is single end full"],
            None,  # This evaluates to meaning running all components
            (487, 7),
            (21, 5),
            14231.17659137577,
            2,  # test data does not have parity between DGE counts and RSEM counts / Encounters rare false positive for check_aggregate_star_unnormalized_counts_table_values_against_samplewise_tables that should only occur when gene counts are zero for multiple strand assessment types
            id="Run all checks",
        ),
    ],
)
def test_updated_protocol_model_single_end(
    glds48_test_dir,
    components,
    data_asset_keys,
    expected_flag_table_shape,
    expected_outlier_table_shape,
    expected_flag_table_fingerprint,
    expected_halt_flag_count,
):
    datasystem = load_data(
        key_sets=data_asset_keys,
        config=("bulkRNASeq", "Latest"),
        root_path=(glds48_test_dir),
        runsheet_path=(glds48_test_dir / "Metadata/GLDS-48_bulkRNASeq_v1_runsheet.csv"),
    )

    report = validate_bulkRNASeq(
        datasystem.dataset,
        report_args={"include_skipped": False},
        protocol_args={"run_components": components},
    )
    assert (
        sum(report["flag_table"]["code_level"] >= 80) == expected_halt_flag_count
    ), f"Found more than expected HALT+ flags: {report['flag_table'].loc[report['flag_table']['code_level'] >= 80].to_dict()}"

    assert (
        report["flag_table"].shape,
        report["outliers"].shape,
        pseudo_fingerprint(report["flag_table"]),
    ) == (
        expected_flag_table_shape,
        expected_outlier_table_shape,
        expected_flag_table_fingerprint,
    )


def test_updated_protocol_model_skipping(glds48_dataSystem):
    report = validate_bulkRNASeq(
        glds48_dataSystem.dataset,
        report_args={"include_skipped": False},
        protocol_args={
            "run_components": [
                "Raw Reads By Sample",
                # "Trimmed Reads By Sample",
                # "STAR Alignments By Sample",
                "Metadata",
                "Raw Reads",
                # "Trim Reads",
            ]
        },
    )

    assert report["flag_table"].shape == (353, 7)
    assert report["outliers"].shape == (1, 1)
    assert pseudo_fingerprint(report["flag_table"]) == 9621.198300283286

    # NOW INCLUDING SKIPPED FLAG TABLE ENTRIES
    # SHOULD MATCH, running all components and not including skips
    report = validate_bulkRNASeq(
        glds48_dataSystem.dataset,
        report_args={"include_skipped": True},
        protocol_args={
            "run_components": [
                "Raw Reads By Sample",
                "Metadata",
                "Raw Reads",
                # "Trim Reads",
            ]
        },
    )

    assert report["flag_table"].shape == (559, 7)
    assert report["outliers"].shape == (1, 1)
    assert pseudo_fingerprint(report["flag_table"]) == 11262.12343470483


def test_updated_protcol_model_printouts_single(glds48_dataSystem):
    vp = validate_bulkRNASeq(glds48_dataSystem.dataset, defer_run=True)

    print(vp.queued_checks(include_individual_checks=False))
    print(vp.queued_checks())


def test_updated_protcol_model_printouts_paired(glds194_dataSystem):
    vp = validate_bulkRNASeq(glds194_dataSystem.dataset, defer_run=True)

    print(vp.queued_checks(include_individual_checks=False))
    print(vp.queued_checks())


def test_report_modification_add_sample_column(glds48_dataSystem):
    report = validate_bulkRNASeq(
        glds48_dataSystem.dataset,
        report_args={"include_skipped": False},
        protocol_args={
            "run_components": [
                "Raw Reads By Sample",
                # "Trimmed Reads By Sample",
                # "STAR Alignments By Sample",
                "Metadata",
                "Raw Reads",
                # "Trim Reads",
            ]
        },
    )
    samples = list(glds48_dataSystem.dataset.samples)
    ValidationProtocol.append_sample_column(report["flag_table"], samples=samples)

    assert report["flag_table"].shape == (353, 8)
    assert report["outliers"].shape == (1, 1)
    assert pseudo_fingerprint(report["flag_table"]) == 9974.198300283286

""" Tests for validation report results, relies on test for loaders passing """
from decimal import DivisionByZero
from pathlib import Path
import os

import pytest
from dp_tools.bulkRNASeq.vv_protocols import (
    validate_bulkRNASeq,
)


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
    "components,expected_flag_table_shape,expected_outlier_table_shape,expected_flag_table_fingerprint",
    [
        pytest.param(
            ["Raw Reads By Sample", "Raw Reads"],
            (145, 7),
            (2, 1),
            3945.0689655172414,
            id="Raw Reads Checks Only",
        ),
        pytest.param(
            ["Trimmed Reads By Sample", "Trim Reads"],
            (171, 7),
            (4, 1),
            4657.116959064328,
            id="Trimmed Reads Checks Only",
        ),
        pytest.param(
            ["STAR Alignments", "STAR Alignments By Sample"],
            (131, 7),
            (6, 1),
            3960.053435114504,
            id="STAR Alignments Checks Only",
        ),
        pytest.param(
            ["RSeQC", "RSeQC By Sample"],
            (109, 7),
            (13, 1),
            3003.3669724770643,
            id="RSeQC Checks Only",
        ),
        pytest.param(
            ["RSEM Counts", "Unnormalized Gene Counts"],
            (4, 7),
            (2, 1),
            178.0,
            id="RSEM Counts Checks Only",
        ),
        pytest.param(
            ["DGE Metadata", "DGE Metadata ERCC", "DGE Output", "DGE Output ERCC"],
            (55, 7),
            (0, 0),
            1566.090909090909,
            id="DGE Checks Only",
        ),
    ],
)
def test_updated_protocol_model_paired_end(
    glds194_dataSystem_STAGE04,
    components,
    expected_flag_table_shape,
    expected_outlier_table_shape,
    expected_flag_table_fingerprint,
):
    report = validate_bulkRNASeq(
        glds194_dataSystem_STAGE04.dataset,
        report_args={"include_skipped": False},
        protocol_args={"run_components": components},
    )

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
    "components,expected_flag_table_shape,expected_outlier_table_shape,expected_flag_table_fingerprint",
    [
        pytest.param(
            ["Raw Reads By Sample", "Raw Reads"],
            (71, 7),
            (1, 1),
            1947.1408450704225,
            id="Raw Reads Checks Only",
        ),
        pytest.param(
            ["Trimmed Reads By Sample", "Trim Reads"],
            (85, 7),
            (2, 1),
            2325.1176470588234,
            id="Trimmed Reads Checks Only",
        ),
        pytest.param(
            ["STAR Alignments", "STAR Alignments By Sample"],
            (141, 7),
            (9, 1),
            3837.0709219858154,
            id="STAR Alignments Checks Only",
        ),
        pytest.param(
            ["RSeQC", "RSeQC By Sample"],
            (88, 7),
            (11, 1),
            2426.340909090909,
            id="RSeQC Checks Only",
        ),
        pytest.param(
            ["RSEM Counts", "Unnormalized Gene Counts"],
            (3, 7),
            (1, 1),
            194.33333333333334,
            id="RSEM Counts Checks Only",
        ),
        pytest.param(
            ["DGE Metadata", "DGE Metadata ERCC", "DGE Output", "DGE Output ERCC"],
            (28, 7),
            (0, 0),
            838.1428571428571,
            id="DGE Checks Only",
        ),
    ],
)
def test_updated_protocol_model_single_end(
    glds48_dataSystem_STAGE04,
    components,
    expected_flag_table_shape,
    expected_outlier_table_shape,
    expected_flag_table_fingerprint,
):
    report = validate_bulkRNASeq(
        glds48_dataSystem_STAGE04.dataset,
        report_args={"include_skipped": False},
        protocol_args={"run_components": components},
    )

    assert (
        report["flag_table"].shape,
        report["outliers"].shape,
        pseudo_fingerprint(report["flag_table"]),
    ) == (
        expected_flag_table_shape,
        expected_outlier_table_shape,
        expected_flag_table_fingerprint,
    )


def test_updated_protocol_model_skipping(glds48_dataSystem_STAGE00):
    report = validate_bulkRNASeq(
        glds48_dataSystem_STAGE00.dataset,
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

    assert report["flag_table"].shape == (72, 7)
    assert report["outliers"].shape == (1, 1)
    assert pseudo_fingerprint(report["flag_table"]) == 1974.138888888889

    # NOW INCLUDING SKIPPED FLAG TABLE ENTRIES
    # SHOULD MATCH, running all components and not including skips
    report = validate_bulkRNASeq(
        glds48_dataSystem_STAGE00.dataset,
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

    assert report["flag_table"].shape == (502, 7)
    assert report["outliers"].shape == (1, 1)
    assert pseudo_fingerprint(report["flag_table"]) == 5397.745019920319


def test_updated_protcol_model_printouts_single(
    glds48_dataSystem_STAGE04
):
    vp = validate_bulkRNASeq(
        glds48_dataSystem_STAGE04.dataset,  defer_run=True
    )

    print(vp.queued_checks(include_individual_checks=False))
    print(vp.queued_checks())
    # 1 / 0  # Manually Validated by inspecting print statement


def test_updated_protcol_model_printouts_paired(
    glds194_dataSystem_STAGE04
):
    vp = validate_bulkRNASeq(
        glds194_dataSystem_STAGE04.dataset,  defer_run=True
    )

    print(vp.queued_checks(include_individual_checks=False))
    print(vp.queued_checks())
    # 1 / 0  # Manually Validated by inspecting print statement

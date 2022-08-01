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
    "data_asset_keys,components,expected_flag_table_shape,expected_outlier_table_shape,expected_flag_table_fingerprint",
    [
        pytest.param(
            ["demuliplexed paired end raw data", "qc reports for paired end raw data"],
            ["Metadata", "Raw Reads By Sample", "Raw Reads"],
            (110, 7),
            (4, 1),
            3070.7272727272725,
            id="Raw Reads Checks Only",
        ),
        pytest.param(
            [
                "paired end trimmed reads",
                "qc reports for paired end trimmed reads data",
            ],
            ["Trimmed Reads By Sample", "Trim Reads"],
            (159, 7),
            (9, 1),
            4333.125786163522,
            id="Trimmed Reads Checks Only",
        ),
        pytest.param(
            ["STAR alignments"],
            ["STAR Alignments", "STAR Alignments By Sample"],
            (159, 7),
            (7, 1),
            5107.968553459119,
            id="STAR Alignments Checks Only",
        ),
        pytest.param(
            ["RSeQC output for paired end data"],
            ["RSeQC", "RSeQC By Sample"],
            (91, 7),
            (13, 1),
            2578.098901098901,
            id="RSeQC Checks Only",
        ),
        pytest.param(
            ["RSEM counts"],
            ["RSEM Counts"],
            (44, 7),
            (3, 1),
            1218.2272727272727,
            id="RSEM Counts Checks Only",
        ),
        pytest.param(
            ["RSEM counts", "STAR alignments"],
            ["Unnormalized Gene Counts"],
            (167, 7),
            (0, 0),
            4559.179640718563,
            id="Unnormalized Gene Counts Checks Only",
        ),
        pytest.param(
            ["DGE Output", "ERCC DGE Output", "RSEM Output"],
            ["DGE Metadata", "DGE Metadata ERCC", "DGE Output", "DGE Output ERCC"],
            (72, 7),
            (0, 0),
            2024.8333333333333,
            id="DGE Checks Only",
        ),
        pytest.param(
            ["is paired end full", "ERCC DGE Output"],
            None,  # This evaluates to meaning running all components
            (650, 7),
            (34, 5),
            18661.67692307692,
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


def test_updated_protcol_model_printouts_single(glds48_dataSystem_STAGE04):
    vp = validate_bulkRNASeq(glds48_dataSystem_STAGE04.dataset, defer_run=True)

    print(vp.queued_checks(include_individual_checks=False))
    print(vp.queued_checks())
    # 1 / 0  # Manually Validated by inspecting print statement


def test_updated_protcol_model_printouts_paired(glds194_dataSystem_STAGE04):
    vp = validate_bulkRNASeq(glds194_dataSystem_STAGE04.dataset, defer_run=True)

    print(vp.queued_checks(include_individual_checks=False))
    print(vp.queued_checks())
    # 1 / 0  # Manually Validated by inspecting print statement


def test_report_modification_add_sample_column(glds48_dataSystem_STAGE00):
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
    samples = [s for s in glds48_dataSystem_STAGE00.dataset.samples]
    ValidationProtocol.append_sample_column(report["flag_table"], samples=samples)

    assert report["flag_table"].shape == (72, 8)
    assert report["outliers"].shape == (1, 1)
    assert pseudo_fingerprint(report["flag_table"]) == 2046.138888888889

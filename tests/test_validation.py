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
            (124, 7),
            (4, 1),
            3478.8870967741937,
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
            id="Trimmed Reads Checks Only",
        ),
        pytest.param(
            ["STAR alignments"],
            ["STAR Alignments", "STAR Alignments By Sample"],
            (175, 7),
            (7, 1),
            5569.685714285714,
            id="STAR Alignments Checks Only",
        ),
        pytest.param(
            ["RSeQC output for paired end data"],
            ["RSeQC", "RSeQC By Sample"],
            (92, 7),
            (13, 1),
            2635.413043478261,
            id="RSeQC Checks Only",
        ),
        pytest.param(
            ["RSEM counts"],
            ["RSEM Counts"],
            (45, 7),
            (3, 1),
            1275.888888888889,
            id="RSEM Counts Checks Only",
        ),
        pytest.param(
            ["RSEM counts", "STAR alignments"],
            ["Unnormalized Gene Counts"],
            (168, 7),
            (0, 0),
            4616.357142857143,
            id="Unnormalized Gene Counts Checks Only",
        ),
        pytest.param(
            ["DGE Output", "ERCC DGE Output", "RSEM Output"],
            ["DGE Metadata", "DGE Metadata ERCC", "DGE Output", "DGE Output ERCC"],
            (73, 7),
            (0, 0),
            2082.232876712329,
            id="DGE Checks Only",
        ),
        pytest.param(
            ["is paired end full", "ERCC DGE Output"],
            None,  # This evaluates to meaning running all components
            (679, 7),
            (34, 5),
            19474.649484536083,
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
    "data_asset_keys,components,expected_flag_table_shape,expected_outlier_table_shape,expected_flag_table_fingerprint",
    [
        pytest.param(
            ["demuliplexed single end raw data", "qc reports for single end raw data"],
            ["Metadata", "Raw Reads By Sample", "Raw Reads"],
            (61, 7),
            (1, 1),
            1738.1475409836066,
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
            id="Trimmed Reads Checks Only",
        ),
        pytest.param(
            ["STAR alignments"],
            ["STAR Alignments", "STAR Alignments By Sample"],
            (187, 7),
            (7, 1),
            5923.545454545455,
            id="STAR Alignments Checks Only",
        ),
        pytest.param(
            ["RSeQC output for single end data"],
            ["RSeQC", "RSeQC By Sample"],
            (66, 7),
            (12, 1),
            1862.909090909091,
            id="RSeQC Checks Only",
        ),
        pytest.param(
            ["RSEM counts"],
            ["RSEM Counts"],
            (47, 7),
            (1, 1),
            1299.2127659574467,
            id="RSEM Counts Checks Only",
        ),
        pytest.param(
            ["RSEM counts", "STAR alignments"],
            ["Unnormalized Gene Counts"],
            (178, 7),
            (0, 0),
            4886.337078651685,
            id="Unnormalized Gene Counts Checks Only",
        ),
        pytest.param(
            ["DGE Output", "RSEM Output"],
            ["DGE Metadata", "DGE Metadata ERCC", "DGE Output", "DGE Output ERCC"],
            (39, 7),
            (0, 0),
            1134.5384615384614,
            id="DGE Checks Only",
        ),
        pytest.param(
            ["is single end full"],
            None,  # This evaluates to meaning running all components
            (487, 7),
            (21, 5),
            14291.29979466119,
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

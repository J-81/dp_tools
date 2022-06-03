""" Tests for validation report results, relies on test for loaders passing """
from collections import defaultdict
from contextlib import contextmanager
from decimal import DivisionByZero
from pathlib import Path
import os
from typing import Callable, TypedDict
from dp_tools.core.check_model import FlagCode, ValidationProtocol

from pytest import MonkeyPatch
import pytest
from dp_tools.bulkRNASeq.entity import BulkRNASeqSample

from dp_tools.bulkRNASeq.loaders import (
    load_BulkRNASeq_STAGE_00,
    load_BulkRNASeq_STAGE_01,
)
from dp_tools.bulkRNASeq.vv_protocols import (
    STAGE,
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


@pytest.fixture
def check_config():
    return Path(os.environ["CONFIG_CHECK_PATH"])


# TODO: Part of an alternative framework not fully implemented

import logging

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)


def test_updated_protocol_model_paired_end_part_1(
    glds194_dataSystem_STAGE04, check_config
):
    report = validate_bulkRNASeq(
        glds194_dataSystem_STAGE04.dataset,
        config_path=check_config,
        report_args={"include_skipped": False},
        protocol_args={
            "run_components": [
                "Raw Reads By Sample",
                "Trimmed Reads By Sample",
                "STAR Alignments By Sample",
                "Metadata",
                "Raw Reads",
                "Trim Reads",
            ]
        },
    )

    assert report["flag_table"].shape == (447, 6)
    assert report["outliers"].shape == (4, 2)
    assert pseudo_fingerprint(report["flag_table"]) == 12463.834451901566


def test_updated_protocol_model_paired_end_custom(
    glds194_dataSystem_STAGE04, glds48_dataSystem_STAGE04, check_config
):
    custom_run_components = {
        "run_components": [
            "DGE Output",
        ]
    }

    """
    report = validate_bulkRNASeq(
        glds194_dataSystem_STAGE04.dataset,
        config_path=check_config,
        report_args={"include_skipped": False},
        protocol_args=custom_run_components,
    )
    """

    print(1)

    report = validate_bulkRNASeq(
        glds48_dataSystem_STAGE04.dataset,
        config_path=check_config,
        report_args={"include_skipped": False},
        protocol_args=custom_run_components,
    )

    print(1)


def test_updated_protocol_model_paired_end_part_2(
    glds194_dataSystem_STAGE04, check_config
):
    report = validate_bulkRNASeq(
        glds194_dataSystem_STAGE04.dataset,
        config_path=check_config,
        report_args={"include_skipped": False},
        protocol_args={
            "run_components": [
                "RSeQC By Sample",
            ]
        },
    )

    assert report["flag_table"].shape == (26, 6)
    # assert report["outliers"].shape == (4, 2)
    # assert pseudo_fingerprint(report["flag_table"]) == 12493.901565995526


def test_updated_protocol_model_single_end(glds48_dataSystem_STAGE04, check_config):
    report = validate_bulkRNASeq(
        glds48_dataSystem_STAGE04.dataset,
        config_path=check_config,
        report_args={"include_skipped": True},
    )

    assert report["flag_table"].shape == (500, 6)
    # now check without skipped entries (should only be the fastq read parity check)
    assert report["flag_table"].loc[
        report["flag_table"].code.apply(lambda v: v.name != "SKIPPED")
    ].shape == (418, 6)
    assert report["outliers"].shape == (13, 5)
    assert pseudo_fingerprint(report["flag_table"]) == 11769.504


def test_updated_protocol_model_paired_end(glds194_dataSystem_STAGE04, check_config):
    report = validate_bulkRNASeq(
        glds194_dataSystem_STAGE04.dataset,
        config_path=check_config,
        report_args={"include_skipped": True},
    )

    assert report["flag_table"].shape == (614, 6)
    # now check without skipped entries (should only be the fastq read parity check)
    assert report["flag_table"].loc[
        report["flag_table"].code.apply(lambda v: v.name != "SKIPPED")
    ].shape == (418, 6)
    assert report["outliers"].shape == (13, 5)
    assert pseudo_fingerprint(report["flag_table"]) == 11769.504


def test_updated_protocol_model_skipping(glds48_dataSystem_STAGE00, check_config):
    report = validate_bulkRNASeq(
        glds48_dataSystem_STAGE00.dataset,
        config_path=check_config,
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

    assert report["flag_table"].shape == (72, 6)
    assert report["outliers"].shape == (1, 1)
    assert pseudo_fingerprint(report["flag_table"]) == 1903.1527777777778

    # NOW INCLUDING SKIPPED FLAG TABLE ENTRIES
    # SHOULD MATCH, running all components and not including skips
    report = validate_bulkRNASeq(
        glds48_dataSystem_STAGE00.dataset,
        config_path=check_config,
        report_args={"include_skipped": True},
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

    assert report["flag_table"].shape == (325, 6)
    assert report["outliers"].shape == (1, 1)
    assert pseudo_fingerprint(report["flag_table"]) == 3659.243076923077


def test_updated_protcol_model_printouts_single(
    glds48_dataSystem_STAGE04, check_config
):
    vp = validate_bulkRNASeq(
        glds48_dataSystem_STAGE04.dataset, config_path=check_config, defer_run=True
    )

    print(vp.queued_checks(include_individual_checks=False))
    print(vp.queued_checks())
    1 / 0  # Manually Validated by inspecting print statement


def test_updated_protcol_model_printouts_paired(
    glds194_dataSystem_STAGE04, check_config
):
    vp = validate_bulkRNASeq(
        glds194_dataSystem_STAGE04.dataset, config_path=check_config, defer_run=True
    )

    print(vp.queued_checks(include_individual_checks=False))
    print(vp.queued_checks())
    1 / 0  # Manually Validated by inspecting print statement

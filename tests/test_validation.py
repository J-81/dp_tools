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


def test_updated_protocol_model_paired_end(glds194_dataSystem_STAGE04, check_config):
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
    assert report["outliers"].shape == (13, 2)
    assert pseudo_fingerprint(report["flag_table"]) == 12574.080536912752


def test_updated_protocol_model_single_end(glds48_dataSystem_STAGE04, check_config):
    report = validate_bulkRNASeq(
        glds48_dataSystem_STAGE04.dataset,
        config_path=check_config,
        report_args={"include_skipped": True},
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

    assert report["flag_table"].shape == (325, 6)
    # now check without skipped entries (should only be the fastq read parity check)
    assert report["flag_table"].loc[
        report["flag_table"].code.apply(lambda v: v.name != "SKIPPED")
    ].shape == (297, 6)
    assert report["outliers"].shape == (14, 2)
    assert pseudo_fingerprint(report["flag_table"]) == 8869.224615384615


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
    assert report["outliers"].shape == (14, 1)
    assert pseudo_fingerprint(report["flag_table"]) == 3679.3046153846153


def test_updated_protcol_model(glds194_dataSystem_STAGE04):
    vp = validate_bulkRNASeq(
        glds194_dataSystem_STAGE04.dataset,
        config_path=check_config,
        defer_run=True,
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

    print(vp.queued_checks())
    # 1/0 Manually Validated by inspecting print statement

""" Tests for validation report results, relies on test for loaders passing """
from decimal import DivisionByZero
from pathlib import Path
import os

from pytest import MonkeyPatch
import pytest
from dp_tools.bulkRNASeq.entity import BulkRNASeqSample

from dp_tools.bulkRNASeq.loaders import (
    load_BulkRNASeq_STAGE_00,
    load_BulkRNASeq_STAGE_01,
)
from dp_tools.bulkRNASeq.vv_protocols import STAGE, BulkRNASeq_VVProtocol


@pytest.fixture(autouse=True)
def mock_dev_exceptions(monkeypatch):
    monkeypatch.setattr(
        "dp_tools.core.check_model.ALLOWED_DEV_EXCEPTIONS", (DivisionByZero)
    )  # ensure unhandled developer exceptions are raised


def test_bulkRNASeq_STAGE00_validation_paired(caplog, glds194_dataSystem_STAGE00):
    """This tests validation as it would be run on dataset after demultiplexing"""
    CAPLEVEL = 20
    caplog.set_level(CAPLEVEL)
    ds = glds194_dataSystem_STAGE00
    vv_protocol = BulkRNASeq_VVProtocol(
        dataset=ds.dataset, dry_run=True, protocol_name="only raw"
    )

    with caplog.at_level(CAPLEVEL):
        vv_protocol.validate_all()

    assert isinstance(vv_protocol.flags["dataset"], dict)
    assert isinstance(vv_protocol.flags["sample"], dict)
    assert isinstance(vv_protocol.flags["component"], dict)

    # second, run with full validation
    with caplog.at_level(CAPLEVEL):
        caplog.clear()
        with MonkeyPatch.context() as m:
            vv_protocol.validate_all()
            df = vv_protocol.flags_to_df()

            df_verbose = vv_protocol.flags_to_df(schema="verbose")

            # assert that no failing flags were raised
            # assert df["flag_code"].max() == 20 # not needed as this tests the truncated data rather than the logic

            # check if appropriate number of flags are raised
            # Currently:
            #   Dataset check : 2
            #   Sample check : 1 per sample
            #   Component checks :
            #       Reads : 1 per component
            assert len(df) == 41
            assert [0] == list(
                df["flag_code"].unique()
            )  # only the dry run code should be returned


def test_bulkRNASeq_STAGE00_validation_paired_no_dry_run(
    caplog, glds194_dataSystem_STAGE00
):
    """This tests validation as it would be run on dataset after demultiplexing"""
    CAPLEVEL = 20
    caplog.set_level(CAPLEVEL)
    ds = glds194_dataSystem_STAGE00
    vv_protocol = BulkRNASeq_VVProtocol(dataset=ds.dataset, protocol_name="only raw")

    with caplog.at_level(CAPLEVEL):
        vv_protocol.validate_all()

    assert isinstance(vv_protocol.flags["dataset"], dict)
    assert isinstance(vv_protocol.flags["sample"], dict)
    assert isinstance(vv_protocol.flags["component"], dict)

    # second, run with full validation
    with caplog.at_level(CAPLEVEL):
        caplog.clear()
        with MonkeyPatch.context() as m:
            vv_protocol.validate_all()
            df = vv_protocol.flags_to_df()

            df_verbose = vv_protocol.flags_to_df(schema="verbose")

            # assert that no failing flags were raised
            # assert df["flag_code"].max() == 20 # not needed as this tests the truncated data rather than the logic

            # check if appropriate number of flags are raised
            # Currently:
            #   Dataset check : 2
            #   Sample check : 1 per sample
            #   Component checks :
            #       Reads : 1 per component
            assert len(df) == 41


def test_bulkRNASeq_STAGE00_validation_paired_with_skips(
    caplog, glds194_dataSystem_STAGE00
):
    """This tests validation as it would be run on dataset after demultiplexing"""
    CAPLEVEL = 20
    caplog.set_level(CAPLEVEL)
    ds = glds194_dataSystem_STAGE00
    vv_protocol = BulkRNASeq_VVProtocol(
        dataset=ds.dataset,
        protocol_name="only raw",
        dry_run=True,
        skip_these_checks={"DATASET_RAWREADS_0001"},
    )

    with caplog.at_level(CAPLEVEL):
        vv_protocol.validate_all()

    assert isinstance(vv_protocol.flags["dataset"], dict)
    assert isinstance(vv_protocol.flags["sample"], dict)
    assert isinstance(vv_protocol.flags["component"], dict)

    # second, run with full validation
    with caplog.at_level(CAPLEVEL):
        caplog.clear()
        with MonkeyPatch.context() as m:
            vv_protocol.validate_all()
            df = vv_protocol.flags_to_df()

            df_verbose = vv_protocol.flags_to_df(schema="verbose")

            # assert that no failing flags were raised
            # assert df["flag_code"].max() == 20 # not needed as this tests the truncated data rather than the logic

            # check if appropriate number of flags are raised
            # Currently:
            #   Dataset check : 2
            #   Sample check : 1 per sample
            #   Component checks :
            #       Reads : 1 per component
            assert len(df) == 41
            assert 0 in df["flag_code"].values  # ensure dry run flag codes returned
            assert 1 in df["flag_code"].values  # ensure skip flag codes returned


def test_bulkRNASeq_STAGE00_validation_paired_with_config(
    caplog, glds194_dataSystem_STAGE00
):
    """This tests validation as it would be run on dataset after demultiplexing"""
    CAPLEVEL = 20
    caplog.set_level(CAPLEVEL)
    ds = glds194_dataSystem_STAGE00
    vv_protocol = BulkRNASeq_VVProtocol(
        dataset=ds.dataset, config=("bulkRNASeq", "0"), protocol_name="only raw"
    )

    with caplog.at_level(CAPLEVEL):
        vv_protocol.validate_all()

    assert isinstance(vv_protocol.flags["dataset"], dict)
    assert isinstance(vv_protocol.flags["sample"], dict)
    assert isinstance(vv_protocol.flags["component"], dict)

    # second, run with full validation
    with caplog.at_level(CAPLEVEL):
        caplog.clear()
        with MonkeyPatch.context() as m:
            vv_protocol.validate_all()
            df = vv_protocol.flags_to_df()

            df_verbose = vv_protocol.flags_to_df(schema="verbose")

            # assert that no failing flags were raised
            # assert df["flag_code"].max() == 20 # not needed as this tests the truncated data rather than the logic

            # check if appropriate number of flags are raised
            # Currently:
            #   Dataset check : 2
            #   Sample check : 1 per sample
            #   Component checks :
            #       Reads : 1 per component
            assert len(df) == 41
            assert 0 not in df["flag_code"].values  # ensure dry run flag codes returned
            assert 1 not in df["flag_code"].values  # ensure skip flag codes returned


def test_bulkRNASeq_STAGE00_validation_single(caplog, glds48_dataSystem_STAGE00):
    """This tests validation as it would be run on dataset after demultiplexing"""
    CAPLEVEL = 20

    caplog.set_level(CAPLEVEL)
    ds = glds48_dataSystem_STAGE00
    vv_protocol = BulkRNASeq_VVProtocol(
        dataset=ds.dataset, protocol_name="only raw", dry_run=True
    )

    with MonkeyPatch.context() as m:
        vv_protocol.validate_all()
        df = vv_protocol.flags_to_df()

        df_verbose = vv_protocol.flags_to_df(schema="verbose")

        # check if appropriate number of flags are raised
        # Currently:
        #   Dataset check : 2
        #   Sample check : 1 per sample
        #   Component checks
        #       Reads : 1 per component (1 per sample)
        assert len(df) == 30


"""
def test_bulkRNASeq_STAGE01_validation_paired(glds194_dataSystem_STAGE01):
    ds = glds194_dataSystem_STAGE01
    vv_protocol = BulkRNASeq_VVProtocol(
        dataset=ds.dataset, stage_names=STAGE.Reads_PreProcessed, dry_run=True
    )

    vv_protocol.validate_all()
    df = vv_protocol.flags_to_df()

    assert len(df) == 81



def test_bulkRNASeq_STAGE01_validation_single(glds48_dataSystem_STAGE01):
    ds = glds48_dataSystem_STAGE01
    vv_protocol = BulkRNASeq_VVProtocol(
        dataset=ds.dataset, stage_names=STAGE.Reads_PreProcessed, dry_run=True
    )

    vv_protocol.validate_all()
    df = vv_protocol.flags_to_df()

    assert len(df) == 59


def test_bulkRNASeq_STAGE02_validation_paired(glds194_dataSystem_STAGE02):
    ds = glds194_dataSystem_STAGE02
    vv_protocol = BulkRNASeq_VVProtocol(
        dataset=ds.dataset, stage_names=STAGE.GenomeAligned, dry_run=True
    )

    vv_protocol.validate_all()
    df = vv_protocol.flags_to_df()

    assert len(df) == 95


def test_bulkRNASeq_STAGE02_validation_single(glds48_dataSystem_STAGE02):
    ds = glds48_dataSystem_STAGE02
    vv_protocol = BulkRNASeq_VVProtocol(
        dataset=ds.dataset, stage_names=STAGE.GenomeAligned, dry_run=True
    )

    vv_protocol.validate_all()
    df = vv_protocol.flags_to_df()

    assert len(df) == 74


def test_bulkRNASeq_STAGE03_validation_paired(glds194_dataSystem_STAGE03):
    ds = glds194_dataSystem_STAGE03
    vv_protocol = BulkRNASeq_VVProtocol(
        dataset=ds.dataset, stage_names=STAGE.GeneCounted, dry_run=True
    )

    vv_protocol.validate_all()
    df = vv_protocol.flags_to_df()

    assert len(df) == 97


def test_bulkRNASeq_STAGE03_validation_single(glds48_dataSystem_STAGE03):
    ds = glds48_dataSystem_STAGE03
    vv_protocol = BulkRNASeq_VVProtocol(
        dataset=ds.dataset, stage_names=STAGE.GeneCounted, dry_run=True
    )

    vv_protocol.validate_all()
    df = vv_protocol.flags_to_df()

    assert len(df) == 76
"""  # DISABLED PENDING REWORK OF ARGS + CONFIG APPROACH


def test_bulkRNASeq_STAGE04_validation_paired(glds194_dataSystem_STAGE04):
    ds = glds194_dataSystem_STAGE04
    vv_protocol = BulkRNASeq_VVProtocol(
        dataset=ds.dataset, protocol_name="full", dry_run=True
    )

    vv_protocol.validate_all()
    df = vv_protocol.flags_to_df()

    assert len(df) == 98


def test_bulkRNASeq_STAGE04_validation_single(glds48_dataSystem_STAGE04):
    ds = glds48_dataSystem_STAGE04
    vv_protocol = BulkRNASeq_VVProtocol(
        dataset=ds.dataset, protocol_name="full", dry_run=True
    )

    vv_protocol.validate_all()
    df = vv_protocol.flags_to_df()

    assert len(df) == 77

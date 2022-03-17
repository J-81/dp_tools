from pathlib import Path
import os

from dp_tools.bulkRNASeq.loaders import (
    load_BulkRNASeq_STAGE_00,
    load_BulkRNASeq_STAGE_01,
)
import pytest

# set for testing
TEST_DIR = Path(os.environ["TEST_ASSETS_DIR"])


@pytest.fixture
def glds194_sample_names():
    return [
        "Mmus_BAL-TAL_LRTN_BSL_Rep1_B7",
        "Mmus_BAL-TAL_RRTN_BSL_Rep2_B8",
        "Mmus_BAL-TAL_RRTN_BSL_Rep3_B9",
        "Mmus_BAL-TAL_RRTN_BSL_Rep4_B10",
        "Mmus_BAL-TAL_LRTN_GC_Rep1_G6",
        "Mmus_BAL-TAL_LRTN_GC_Rep2_G8",
        "Mmus_BAL-TAL_LRTN_GC_Rep3_G9",
        "Mmus_BAL-TAL_RRTN_GC_Rep4_G10",
        "Mmus_BAL-TAL_LRTN_FLT_Rep1_F6",
        "Mmus_BAL-TAL_LRTN_FLT_Rep2_F7",
        "Mmus_BAL-TAL_LRTN_FLT_Rep3_F8",
        "Mmus_BAL-TAL_LRTN_FLT_Rep4_F9",
        "Mmus_BAL-TAL_LRTN_FLT_Rep5_F10",
    ]


@pytest.fixture
def glds48_sample_names():
    return [
        "Mmus_C57-6J_LVR_GC_I_Rep1_M31",
        "Mmus_C57-6J_LVR_GC_I_Rep2_M32",
        "Mmus_C57-6J_LVR_FLT_I_Rep1_M21",
        "Mmus_C57-6J_LVR_FLT_I_Rep2_M22",
        "Mmus_C57-6J_LVR_GC_C_Rep1_M36",
        "Mmus_C57-6J_LVR_GC_C_Rep2_M37",
        "Mmus_C57-6J_LVR_GC_C_Rep3_M38",
        "Mmus_C57-6J_LVR_GC_C_Rep4_M39",
        "Mmus_C57-6J_LVR_GC_C_Rep5_M40",
        "Mmus_C57-6J_LVR_FLT_C_Rep1_M25",
        "Mmus_C57-6J_LVR_FLT_C_Rep2_M26",
        "Mmus_C57-6J_LVR_FLT_C_Rep3_M27",
        "Mmus_C57-6J_LVR_FLT_C_Rep4_M28",
        "Mmus_C57-6J_LVR_FLT_C_Rep5_M30",
    ]


@pytest.fixture
def glds194_test_dir():
    return TEST_DIR / "GLDS-194_TruncatedProcessed"


@pytest.fixture
def glds48_test_dir():
    return TEST_DIR / "GLDS-48_TruncatedProcessed"


@pytest.fixture
def typo_test_dir():
    return TEST_DIR / "GLDS-48_BUTWITHTYPOS"


def test_bulkRNASeq_with_bad_rootdir(typo_test_dir):

    with pytest.raises(FileNotFoundError):
        ds = load_BulkRNASeq_STAGE_00(typo_test_dir, dataSystem_name="GLDS-48")


def test_bulkRNASeq_STAGE00_single(caplog, glds48_test_dir, glds48_sample_names):
    """ Tests loader for state after demultiplexing for single end study """

    caplog.set_level(0)
    ds = load_BulkRNASeq_STAGE_00(glds48_test_dir, dataSystem_name="GLDS-48")

    # pull dataset
    dataset = ds.datasets["GLDS-48__BulkRNASeq"]

    assert list(dataset.samples.keys()) == glds48_sample_names

    # check expected loaded components [raw directory]
    for sample in ds.dataset.samples.values():
        assert len(list(sample.rawReads.multiQCDir.path.iterdir())) == 2


def test_bulkRNASeq_STAGE00_paired(caplog, glds194_test_dir, glds194_sample_names):
    """ Tests loader for state after demultiplexing for single end study """

    caplog.set_level(0)
    ds = load_BulkRNASeq_STAGE_00(glds194_test_dir, dataSystem_name="GLDS-194")

    # pull dataset
    dataset = ds.datasets["GLDS-194__BulkRNASeq"]

    assert list(dataset.samples) == glds194_sample_names

    # check expected loaded components [raw directory]
    for sample in ds.dataset.samples.values():
        assert len(list(sample.rawForwardReads.multiQCDir.path.iterdir())) == 2
        assert len(list(sample.rawReverseReads.multiQCDir.path.iterdir())) == 2


def test_bulkRNASeq_STAGE01_single(caplog, glds48_test_dir, glds48_sample_names):
    """ Tests loader for state after demultiplexing for single end study """

    caplog.set_level(0)
    ds = load_BulkRNASeq_STAGE_01(
        *load_BulkRNASeq_STAGE_00(
            glds48_test_dir, dataSystem_name="GLDS-48", stack=True
        )
    )

    # pull dataset
    dataset = ds.datasets["GLDS-48__BulkRNASeq"]

    assert list(dataset.samples.keys()) == glds48_sample_names

    # check expected loaded components [raw directory]
    for sample in ds.dataset.samples.values():
        assert len(list(sample.rawReads.multiQCDir.path.iterdir())) == 2
        assert len(list(sample.trimReads.multiQCDir.path.iterdir())) == 2


def test_bulkRNASeq_STAGE01_paired(caplog, glds194_test_dir):
    """ Tests loader for state after demultiplexing for single end study """

    caplog.set_level(0)
    ds = load_BulkRNASeq_STAGE_01(
        *load_BulkRNASeq_STAGE_00(
            glds194_test_dir, dataSystem_name="GLDS-194", stack=True
        )
    )

    # pull dataset
    dataset = ds.datasets["GLDS-194__BulkRNASeq"]

    assert list(dataset.samples.keys()) == [
        "Mmus_BAL-TAL_LRTN_BSL_Rep1_B7",
        "Mmus_BAL-TAL_RRTN_BSL_Rep2_B8",
        "Mmus_BAL-TAL_RRTN_BSL_Rep3_B9",
        "Mmus_BAL-TAL_RRTN_BSL_Rep4_B10",
        "Mmus_BAL-TAL_LRTN_GC_Rep1_G6",
        "Mmus_BAL-TAL_LRTN_GC_Rep2_G8",
        "Mmus_BAL-TAL_LRTN_GC_Rep3_G9",
        "Mmus_BAL-TAL_RRTN_GC_Rep4_G10",
        "Mmus_BAL-TAL_LRTN_FLT_Rep1_F6",
        "Mmus_BAL-TAL_LRTN_FLT_Rep2_F7",
        "Mmus_BAL-TAL_LRTN_FLT_Rep3_F8",
        "Mmus_BAL-TAL_LRTN_FLT_Rep4_F9",
        "Mmus_BAL-TAL_LRTN_FLT_Rep5_F10",
    ]

    # check expected loaded components [raw directory]
    for sample in ds.dataset.samples.values():
        assert len(list(sample.rawForwardReads.multiQCDir.path.iterdir())) == 2
        assert len(list(sample.rawReverseReads.multiQCDir.path.iterdir())) == 2
        assert len(list(sample.trimForwardReads.multiQCDir.path.iterdir())) == 2
        assert len(list(sample.trimReverseReads.multiQCDir.path.iterdir())) == 2


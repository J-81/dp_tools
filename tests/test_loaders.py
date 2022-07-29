import logging
from pathlib import Path

from dp_tools.bulkRNASeq.loaders import (
    load_BulkRNASeq_STAGE_00,
    load_BulkRNASeq_STAGE_01,
    load_BulkRNASeq_STAGE_02,
    load_BulkRNASeq_STAGE_0201,
    load_BulkRNASeq_STAGE_03,
    load_BulkRNASeq_STAGE_04,
)

from dp_tools.core.loaders import load_data
from dp_tools.core.configuration import available_data_asset_keys

import pytest


def test_bulkRNASeq_with_bad_rootdir(typo_test_dir):

    with pytest.raises(FileNotFoundError):
        ds = load_BulkRNASeq_STAGE_00(typo_test_dir, dataSystem_name="GLDS-48")


def test_bulkRNASeq_STAGE00_single(caplog, glds48_test_dir, glds48_sample_names):
    """Tests loader for state after demultiplexing for single end study"""

    caplog.set_level(0)
    ds = load_BulkRNASeq_STAGE_00(
        glds48_test_dir, dataSystem_name="GLDS-48", config="Latest"
    )

    # pull dataset
    dataset = ds.datasets["GLDS-48__BulkRNASeq"]

    assert list(dataset.samples.keys()) == glds48_sample_names


def test_bulkRNASeq_STAGE00_paired(caplog, glds194_test_dir, glds194_sample_names):
    """Tests loader for state after demultiplexing for single end study"""

    caplog.set_level(0)
    ds = load_BulkRNASeq_STAGE_00(
        glds194_test_dir, config="Latest", dataSystem_name="GLDS-194"
    )

    # pull dataset
    dataset = ds.datasets["GLDS-194__BulkRNASeq"]

    assert list(dataset.samples) == glds194_sample_names


def test_bulkRNASeq_STAGE01_paired(caplog, glds194_test_dir, glds194_sample_names):
    """Tests loader for state after demultiplexing for single end study"""

    caplog.set_level(0)
    ds = load_BulkRNASeq_STAGE_01(
        *load_BulkRNASeq_STAGE_00(
            glds194_test_dir, dataSystem_name="GLDS-194", config="Latest", stack=True
        ),
        config="Latest",
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


def test_bulkRNASeq_STAGE01_single(caplog, glds48_test_dir, glds48_sample_names):
    """Tests loader for state after demultiplexing for single end study"""

    caplog.set_level(0)
    ds = load_BulkRNASeq_STAGE_01(
        *load_BulkRNASeq_STAGE_00(
            glds48_test_dir, config="Latest", dataSystem_name="GLDS-48", stack=True
        ),
        config="Latest",
    )

    # pull dataset
    dataset = ds.datasets["GLDS-48__BulkRNASeq"]

    assert list(dataset.samples.keys()) == glds48_sample_names


def test_bulkRNASeq_STAGE02_paired(caplog, glds194_test_dir):
    ds = load_BulkRNASeq_STAGE_02(
        *load_BulkRNASeq_STAGE_01(
            *load_BulkRNASeq_STAGE_00(
                glds194_test_dir,
                config="Latest",
                dataSystem_name="GLDS-194",
                stack=True,
            ),
            config="Latest",
            stack=True,
        ),
        config="Latest",
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


def test_bulkRNASeq_STAGE02_single(caplog, glds48_test_dir, glds48_sample_names):
    ds = load_BulkRNASeq_STAGE_02(
        *load_BulkRNASeq_STAGE_01(
            *load_BulkRNASeq_STAGE_00(
                glds48_test_dir, config="Latest", dataSystem_name="GLDS-48", stack=True
            ),
            config="Latest",
            stack=True,
        ),
        config="Latest",
    )

    # pull dataset
    dataset = ds.datasets["GLDS-48__BulkRNASeq"]

    assert list(dataset.samples.keys()) == glds48_sample_names


def test_bulkRNASeq_STAGE0201_paired(caplog, glds194_test_dir):
    ds = load_BulkRNASeq_STAGE_0201(
        *load_BulkRNASeq_STAGE_02(
            *load_BulkRNASeq_STAGE_01(
                *load_BulkRNASeq_STAGE_00(
                    glds194_test_dir,
                    config="Latest",
                    dataSystem_name="GLDS-194",
                    stack=True,
                ),
                config="Latest",
                stack=True,
            ),
            config="Latest",
            stack=True,
        ),
        config="Latest",
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


def test_bulkRNASeq_STAGE0201_single(caplog, glds48_test_dir, glds48_sample_names):
    ds = load_BulkRNASeq_STAGE_0201(
        *load_BulkRNASeq_STAGE_02(
            *load_BulkRNASeq_STAGE_01(
                *load_BulkRNASeq_STAGE_00(
                    glds48_test_dir,
                    config="Latest",
                    dataSystem_name="GLDS-48",
                    stack=True,
                ),
                config="Latest",
                stack=True,
            ),
            config="Latest",
            stack=True,
        ),
        config="Latest",
    )

    # pull dataset
    dataset = ds.datasets["GLDS-48__BulkRNASeq"]

    assert list(dataset.samples.keys()) == glds48_sample_names


def test_bulkRNASeq_STAGE03_paired(caplog, glds194_test_dir, glds194_sample_names):
    ds = load_BulkRNASeq_STAGE_03(
        *load_BulkRNASeq_STAGE_0201(
            *load_BulkRNASeq_STAGE_02(
                *load_BulkRNASeq_STAGE_01(
                    *load_BulkRNASeq_STAGE_00(
                        glds194_test_dir,
                        config="Latest",
                        dataSystem_name="GLDS-194",
                        stack=True,
                    ),
                    config="Latest",
                    stack=True,
                ),
                config="Latest",
                stack=True,
            ),
            config="Latest",
            stack=True,
        ),
        config="Latest",
    )

    # pull dataset
    dataset = ds.datasets["GLDS-194__BulkRNASeq"]

    assert list(dataset.samples.keys()) == glds194_sample_names


def test_bulkRNASeq_STAGE03_single(caplog, glds48_test_dir, glds48_sample_names):
    ds = load_BulkRNASeq_STAGE_03(
        *load_BulkRNASeq_STAGE_0201(
            *load_BulkRNASeq_STAGE_02(
                *load_BulkRNASeq_STAGE_01(
                    *load_BulkRNASeq_STAGE_00(
                        glds48_test_dir,
                        config="Latest",
                        dataSystem_name="GLDS-48",
                        stack=True,
                    ),
                    config="Latest",
                    stack=True,
                ),
                config="Latest",
                stack=True,
            ),
            config="Latest",
            stack=True,
        ),
        config="Latest",
    )

    # pull dataset
    dataset = ds.datasets["GLDS-48__BulkRNASeq"]

    assert list(dataset.samples.keys()) == glds48_sample_names


def test_bulkRNASeq_STAGE04_paired(caplog, glds194_test_dir, glds194_sample_names):
    ds = load_BulkRNASeq_STAGE_04(
        *load_BulkRNASeq_STAGE_03(
            *load_BulkRNASeq_STAGE_0201(
                *load_BulkRNASeq_STAGE_02(
                    *load_BulkRNASeq_STAGE_01(
                        *load_BulkRNASeq_STAGE_00(
                            glds194_test_dir,
                            config="Latest",
                            dataSystem_name="GLDS-194",
                            stack=True,
                        ),
                        config="Latest",
                        stack=True,
                    ),
                    config="Latest",
                    stack=True,
                ),
                config="Latest",
                stack=True,
            ),
            config="Latest",
            stack=True,
        ),
        config="Latest",
    )

    # pull dataset
    dataset = ds.datasets["GLDS-194__BulkRNASeq"]

    assert list(dataset.samples) == glds194_sample_names


def test_bulkRNASeq_STAGE04_single(caplog, glds48_test_dir, glds48_sample_names):
    ds = load_BulkRNASeq_STAGE_04(
        *load_BulkRNASeq_STAGE_03(
            *load_BulkRNASeq_STAGE_0201(
                *load_BulkRNASeq_STAGE_02(
                    *load_BulkRNASeq_STAGE_01(
                        *load_BulkRNASeq_STAGE_00(
                            glds48_test_dir,
                            config="Latest",
                            dataSystem_name="GLDS-48",
                            stack=True,
                        ),
                        config="Latest",
                        stack=True,
                    ),
                    config="Latest",
                    stack=True,
                ),
                config="Latest",
                stack=True,
            ),
            config="Latest",
            stack=True,
        ),
        config="Latest",
    )

    # pull dataset
    dataset = ds.datasets["GLDS-48__BulkRNASeq"]

    assert list(dataset.samples.keys()) == glds48_sample_names


def test_bulkRNASeq_STAGE04_single_global_no_validation(
    caplog, glds48_test_dir, glds48_sample_names
):
    ds = load_BulkRNASeq_STAGE_04(
        *load_BulkRNASeq_STAGE_03(
            *load_BulkRNASeq_STAGE_0201(
                *load_BulkRNASeq_STAGE_02(
                    *load_BulkRNASeq_STAGE_01(
                        *load_BulkRNASeq_STAGE_00(
                            glds48_test_dir,
                            config="Latest",
                            dataSystem_name="GLDS-48",
                            stack=True,
                            validation_enabled=False,
                        ),
                        config="Latest",
                        stack=True,
                        validation_enabled=False,
                    ),
                    config="Latest",
                    stack=True,
                    validation_enabled=False,
                ),
                config="Latest",
                stack=True,
                validation_enabled=False,
            ),
            config="Latest",
            stack=True,
            validation_enabled=False,
        ),
        config="Latest",
        validation_enabled=False,
    )

    # pull dataset
    dataset = ds.datasets["GLDS-48__BulkRNASeq"]

    assert list(dataset.samples.keys()) == glds48_sample_names


def test_list_keys_in_config():
    keys = available_data_asset_keys(config=("bulkRNASeq", "1"))
    assert len(keys) == 74


def test_key_based_loader_GLDS48(root_test_dir):

    ds = load_data(
        key_sets=["is single end full"],
        config=("bulkRNASeq", "Latest"),
        root_path=(root_test_dir / "GLDS-48"),
        runsheet_path=(
            root_test_dir / "GLDS-48/Metadata/GLDS-48_bulkRNASeq_v1_runsheet.csv"
        ),
    )


def test_key_based_loader_GLDS194(root_test_dir, caplog):
    caplog.set_level(logging.INFO)

    load_data(
        key_sets=["is paired end full", "has ercc"],
        config=("bulkRNASeq", "Latest"),
        root_path=(root_test_dir / "GLDS-194"),
        runsheet_path=(
            root_test_dir / "GLDS-194/Metadata/GLDS-194_bulkRNASeq_v1_runsheet.csv"
        ),
    )


def test_key_based_loader_bad_keys(root_test_dir):
    with pytest.raises(AssertionError):
        load_data(
            keys=["typoassetkey"],
            config=("bulkRNASeq", "Latest"),
            root_path=(root_test_dir / "GLDS-194"),
            runsheet_path=(
                root_test_dir / "GLDS-194/Metadata/GLDS-194_bulkRNASeq_v1_runsheet.csv"
            ),
        )

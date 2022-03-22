""" Tests for specific validation checks """
from unittest.mock import MagicMock
from pytest import MonkeyPatch
import pytest

from dp_tools.bulkRNASeq.checks import (
    DATASET_GENOMEALIGNMENTS_0001,
    DATASET_METADATA_0001,
    DATASET_RAWREADS_0001,
    DATASET_TRIMREADS_0001,
)


@pytest.fixture(autouse=True)
def mock_dev_exceptions(monkeypatch):
    monkeypatch.setattr(
        "dp_tools.core.check_model.ALLOWED_DEV_EXCEPTIONS", (SystemExit)
    )  # ensure unhandled developer exceptions are raised


@pytest.fixture(scope="module")
def test_DATASET_METADATA_0001():
    return DATASET_METADATA_0001()


@pytest.fixture(scope="module")
def test_DATASET_RAWREADS_0001():
    return DATASET_RAWREADS_0001()


@pytest.fixture(scope="module")
def test_DATASET_TRIMREADS_0001():
    return DATASET_TRIMREADS_0001()


@pytest.fixture(scope="module")
def test_DATASET_GENOMEALIGNMENTS_0001():
    return DATASET_GENOMEALIGNMENTS_0001()


def test_DATASET_METADATA_0001_paired(
    glds194_dataSystem_STAGE00, test_DATASET_METADATA_0001
):
    dataset = glds194_dataSystem_STAGE00.dataset
    testCheck = test_DATASET_METADATA_0001

    # expected GREEN
    flag = testCheck.validate(dataset)
    assert flag.code.name == "GREEN"
    assert (
        flag.message
        == "All expected metadata is accessible and populated. {'paired_end': True, 'has_ercc': True}"
    )

    # expected HALT1
    with MonkeyPatch.context() as m:
        m.setattr(dataset.metadata, "has_ercc", None)
        flag = testCheck.validate(dataset)
        assert flag.code.name == "HALT1"
        assert flag.message == "Missing expected metadata fields: ['has_ercc']"


def test_DATASET_METADATA_0001_single(
    glds48_dataSystem_STAGE00, test_DATASET_METADATA_0001
):
    dataset = glds48_dataSystem_STAGE00.dataset
    testCheck = test_DATASET_METADATA_0001

    # expected GREEN
    flag = testCheck.validate(dataset)
    assert flag.code.name == "GREEN"
    assert (
        flag.message
        == "All expected metadata is accessible and populated. {'paired_end': False, 'has_ercc': False}"
    )

    # expected HALT1
    with MonkeyPatch.context() as m:
        m.setattr(dataset.metadata, "has_ercc", None)
        flag = testCheck.validate(dataset)
        assert flag.code.name == "HALT1"
        assert flag.message == "Missing expected metadata fields: ['has_ercc']"


def test_DATASET_RAWREADS_0001_paired(
    glds194_dataSystem_STAGE00, test_DATASET_RAWREADS_0001
):
    dataset = glds194_dataSystem_STAGE00.dataset
    testCheck = test_DATASET_RAWREADS_0001

    assert (
        testCheck.description
        == "Check that the reads stats (source from FastQC) have no outliers among samples for the following metrics: ['percent_gc', 'avg_sequence_length', 'total_sequences', 'percent_duplicates']. Yellow Flagged Outliers are defined as a being 1 standard deviations away from the mean. Red Flagged Outliers are defined as a being 2 standard deviations away from the mean. "
    )

    # expected YELLOW1
    with MonkeyPatch.context() as m:
        m.setitem(testCheck.config, "red_standard_deviation_threshold", 10)
        flag = testCheck.validate(dataset)
        assert flag.code.name == "YELLOW1"
        print(flag.message)
        assert flag.message.startswith(
            "Outliers detected as follows (values are rounded number of standard deviations from middle):"
        )
        assert flag.message_args["outliers"] == {
            "percent_gc": {
                "Mmus_BAL-TAL_RRTN_BSL_Rep4_B10:rawForwardReads": -2.1527424662749164,
                "Mmus_BAL-TAL_LRTN_GC_Rep1_G6:rawForwardReads": 1.0456177693335293,
                "Mmus_BAL-TAL_LRTN_GC_Rep2_G8:rawForwardReads": 1.8452078282356408,
                "Mmus_BAL-TAL_LRTN_FLT_Rep2_F7:rawForwardReads": 1.0456177693335293,
                "Mmus_BAL-TAL_RRTN_BSL_Rep4_B10:rawReverseReads": -2.591001102525186,
                "Mmus_BAL-TAL_LRTN_GC_Rep2_G8:rawReverseReads": 1.6193756890782367,
            },
            "percent_duplicates": {
                "Mmus_BAL-TAL_RRTN_BSL_Rep2_B8:rawForwardReads": 1.060991783534611,
                "Mmus_BAL-TAL_RRTN_BSL_Rep4_B10:rawForwardReads": -2.3872315129528747,
                "Mmus_BAL-TAL_LRTN_GC_Rep3_G9:rawForwardReads": 1.060991783534611,
                "Mmus_BAL-TAL_RRTN_GC_Rep4_G10:rawForwardReads": 1.060991783534611,
                "Mmus_BAL-TAL_LRTN_FLT_Rep4_F9:rawForwardReads": -1.0079421943578806,
                "Mmus_BAL-TAL_LRTN_FLT_Rep5_F10:rawForwardReads": 1.060991783534611,
                "Mmus_BAL-TAL_RRTN_BSL_Rep4_B10:rawReverseReads": -1.4366010720592737,
                "Mmus_BAL-TAL_RRTN_GC_Rep4_G10:rawReverseReads": 1.4814948555611258,
                "Mmus_BAL-TAL_LRTN_FLT_Rep3_F8:rawReverseReads": 1.4814948555611258,
                "Mmus_BAL-TAL_LRTN_FLT_Rep4_F9:rawReverseReads": -1.4366010720592737,
            },
        }

    # expected RED1
    with MonkeyPatch.context() as m:
        mp_sample = list(dataset.samples.values())[0]
        m.setitem(
            mp_sample.rawForwardReads.mqcData["FastQC"]["General_Stats"],
            "total_sequences",
            9999999999999,
        )
        flag = testCheck.validate(dataset)
        assert flag.code.name == "RED1"
        print(flag.message)
        assert flag.message.startswith(
            "Outliers detected as follows (values are rounded number of standard deviations from middle):"
        )
        assert flag.message_args["outliers"] == {
            "percent_gc": {
                "Mmus_BAL-TAL_RRTN_BSL_Rep4_B10:rawForwardReads": -2.1527424662749164,
                "Mmus_BAL-TAL_LRTN_GC_Rep1_G6:rawForwardReads": 1.0456177693335293,
                "Mmus_BAL-TAL_LRTN_GC_Rep2_G8:rawForwardReads": 1.8452078282356408,
                "Mmus_BAL-TAL_LRTN_FLT_Rep2_F7:rawForwardReads": 1.0456177693335293,
                "Mmus_BAL-TAL_RRTN_BSL_Rep4_B10:rawReverseReads": -2.591001102525186,
                "Mmus_BAL-TAL_LRTN_GC_Rep2_G8:rawReverseReads": 1.6193756890782367,
            },
            "total_sequences": {
                "Mmus_BAL-TAL_LRTN_BSL_Rep1_B7:rawForwardReads": 3.328201177351375
            },
            "percent_duplicates": {
                "Mmus_BAL-TAL_RRTN_BSL_Rep2_B8:rawForwardReads": 1.060991783534611,
                "Mmus_BAL-TAL_RRTN_BSL_Rep4_B10:rawForwardReads": -2.3872315129528747,
                "Mmus_BAL-TAL_LRTN_GC_Rep3_G9:rawForwardReads": 1.060991783534611,
                "Mmus_BAL-TAL_RRTN_GC_Rep4_G10:rawForwardReads": 1.060991783534611,
                "Mmus_BAL-TAL_LRTN_FLT_Rep4_F9:rawForwardReads": -1.0079421943578806,
                "Mmus_BAL-TAL_LRTN_FLT_Rep5_F10:rawForwardReads": 1.060991783534611,
                "Mmus_BAL-TAL_RRTN_BSL_Rep4_B10:rawReverseReads": -1.4366010720592737,
                "Mmus_BAL-TAL_RRTN_GC_Rep4_G10:rawReverseReads": 1.4814948555611258,
                "Mmus_BAL-TAL_LRTN_FLT_Rep3_F8:rawReverseReads": 1.4814948555611258,
                "Mmus_BAL-TAL_LRTN_FLT_Rep4_F9:rawReverseReads": -1.4366010720592737,
            },
        }


def test_DATASET_RAWREADS_0001_single(
    glds48_dataSystem_STAGE00, test_DATASET_RAWREADS_0001
):
    dataset = glds48_dataSystem_STAGE00.dataset
    testCheck = test_DATASET_RAWREADS_0001

    assert (
        testCheck.description
        == "Check that the reads stats (source from FastQC) have no outliers among samples for the following metrics: ['percent_gc', 'avg_sequence_length', 'total_sequences', 'percent_duplicates']. Yellow Flagged Outliers are defined as a being 1 standard deviations away from the mean. Red Flagged Outliers are defined as a being 2 standard deviations away from the mean. "
    )

    # expected YELLOW1
    with MonkeyPatch.context() as m:
        m.setitem(testCheck.config, "red_standard_deviation_threshold", 10)
        flag = testCheck.validate(dataset)
        assert flag.code.name == "YELLOW1"
        print(flag.message)
        assert flag.message.startswith(
            "Outliers detected as follows (values are rounded number of standard deviations from middle):"
        )
        assert flag.message_args["outliers"] == {
            "percent_gc": {
                "Mmus_C57-6J_LVR_GC_I_Rep2_M32": 1.7138966256236114,
                "Mmus_C57-6J_LVR_FLT_I_Rep1_M21": 1.1685658811070083,
                "Mmus_C57-6J_LVR_FLT_I_Rep2_M22": 1.7138966256236114,
                "Mmus_C57-6J_LVR_FLT_C_Rep1_M25": -1.0127570969594042,
                "Mmus_C57-6J_LVR_FLT_C_Rep2_M26": -1.5580878414760073,
                "Mmus_C57-6J_LVR_FLT_C_Rep4_M28": -1.0127570969594042,
            },
            "percent_duplicates": {
                "Mmus_C57-6J_LVR_GC_I_Rep1_M31": -1.2605690005832242,
                "Mmus_C57-6J_LVR_FLT_I_Rep1_M21": -1.2605690005832242,
                "Mmus_C57-6J_LVR_GC_C_Rep1_M36": 1.0271302967715157,
                "Mmus_C57-6J_LVR_FLT_C_Rep2_M26": 1.0271302967715157,
                "Mmus_C57-6J_LVR_FLT_C_Rep4_M28": 1.6807586674442985,
                "Mmus_C57-6J_LVR_FLT_C_Rep5_M30": 1.3539444821079072,
            },
        }

    # expected RED1
    with MonkeyPatch.context() as m:
        mp_sample = list(dataset.samples.values())[0]
        m.setitem(
            mp_sample.rawReads.mqcData["FastQC"]["General_Stats"],
            "total_sequences",
            9999999999999,
        )
        flag = testCheck.validate(dataset)
        assert flag.code.name == "RED1"
        print(flag.message)
        assert flag.message.startswith(
            "Outliers detected as follows (values are rounded number of standard deviations from middle):"
        )
        assert flag.message_args["outliers"] == {'percent_gc': {'Mmus_C57-6J_LVR_GC_I_Rep2_M32': 1.7138966256236114, 'Mmus_C57-6J_LVR_FLT_I_Rep1_M21': 1.1685658811070083, 'Mmus_C57-6J_LVR_FLT_I_Rep2_M22': 1.7138966256236114, 'Mmus_C57-6J_LVR_FLT_C_Rep1_M25': -1.0127570969594042, 'Mmus_C57-6J_LVR_FLT_C_Rep2_M26': -1.5580878414760073, 'Mmus_C57-6J_LVR_FLT_C_Rep4_M28': -1.0127570969594042}, 'total_sequences': {'Mmus_C57-6J_LVR_GC_I_Rep1_M31': 3.4743961448615166}, 'percent_duplicates': {'Mmus_C57-6J_LVR_GC_I_Rep1_M31': -1.2605690005832242, 'Mmus_C57-6J_LVR_FLT_I_Rep1_M21': -1.2605690005832242, 'Mmus_C57-6J_LVR_GC_C_Rep1_M36': 1.0271302967715157, 'Mmus_C57-6J_LVR_FLT_C_Rep2_M26': 1.0271302967715157, 'Mmus_C57-6J_LVR_FLT_C_Rep4_M28': 1.6807586674442985, 'Mmus_C57-6J_LVR_FLT_C_Rep5_M30': 1.3539444821079072}}


def test_DATASET_TRIMREADS_0001_paired():
    raise NotImplementedError


def test_DATASET_TRIMREADS_0001_single():
    raise NotImplementedError


def test_DATASET_GENOMEALIGNMENTS_0001_paired():
    raise NotImplementedError


def test_DATASET_GENOMEALIGNMENTS_0001_single():
    raise NotImplementedError


def test_DATASET_RSEQCANALYSIS_0001_paired():
    raise NotImplementedError


def test_DATASET_RSEQCANALYSIS_0001_single():
    raise NotImplementedError


def test_DATASET_GENECOUNTS_0001_paired():
    raise NotImplementedError


def test_DATASET_GENECOUNTS_0001_single():
    raise NotImplementedError


def test_DATASET_DIFFERENTIALGENEEXPRESSION_0001_paired():
    raise NotImplementedError


def test_DATASET_DIFFERENTIALGENEEXPRESSION_0001_single():
    raise NotImplementedError

""" Tests for specific validation checks """
from unittest.mock import MagicMock
from pytest import MonkeyPatch
import pytest

from dp_tools.bulkRNASeq.checks import DATASET_METADATA_0001


@pytest.fixture(scope="module")
def test_DATASET_METADATA_0001():
    return DATASET_METADATA_0001()


@pytest.fixture(scope="module")
def test_DATASET_RAWREADS_0001():
    return DATASET_RAWREADS_0001()


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
    glds194_dataSystem_STAGE00, test_DATASET_RAWREADS_0001_paired
):
    raise NotImplementedError


def test_DATASET_RAWREADS_0001_single(
    glds194_dataSystem_STAGE00, test_DATASET_RAWREADS_0001_paired
):
    raise NotImplementedError


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

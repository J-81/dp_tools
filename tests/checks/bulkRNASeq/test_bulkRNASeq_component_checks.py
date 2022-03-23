import copy
import gzip
from io import StringIO
import io
import os
from pathlib import Path
from unittest.mock import MagicMock
import pickle

import pytest
from dp_tools.bulkRNASeq.checks import COMPONENT_RAWREADS_0001, COMPONENT_TRIMREADS_0001
from dp_tools.bulkRNASeq.loaders import load_BulkRNASeq_STAGE_00
from dp_tools.core import check_model
from dp_tools.core.entity_model import DataFile

from pytest import MonkeyPatch


@pytest.fixture(scope="module")
def test_COMPONENT_RAWREADS_0001():
    return COMPONENT_RAWREADS_0001()


@pytest.fixture(scope="module")
def test_COMPONENT_TRIMREADS_0001():
    return COMPONENT_TRIMREADS_0001()

@pytest.fixture(autouse=True)
def mock_dev_exceptions(monkeypatch):
    monkeypatch.setattr(
        "dp_tools.core.check_model.ALLOWED_DEV_EXCEPTIONS", (SystemExit)
    )  # ensure unhandled developer exceptions are raised

def test_COMPONENT_RAWREADS_0001_paired(
    glds194_dataSystem_STAGE00, test_COMPONENT_RAWREADS_0001
):
    ds = glds194_dataSystem_STAGE00
    testCheck = test_COMPONENT_RAWREADS_0001
    test_component = list(ds.dataset.samples.values())[0].rawForwardReads

    # expected GREEN
    flag = testCheck.validate(test_component)
    assert flag.maxCode.name == "GREEN"

    # expected HALT1
    with MonkeyPatch.context() as m:
        m.setattr(test_component, "fastqGZ", None)
        m.setattr(test_component, "fastqcReportZIP", None)
        flag = testCheck.validate(test_component)
        assert flag.maxCode.name == "HALT1"
        assert flag.message == "Missing expected files: ['fastqGZ', 'fastqcReportZIP']"

    # expected HALT2
    # TODO: requires a true truncated test file or something similar
    # with MonkeyPatch.context() as m:
    #    m.setattr(gzip, "open", lambda: raise EOFError)
    #    flag = COMPONENT_RAWREADS_0001.validate(test_component)
    #    assert flag.maxCode.name == 'HALT2'

    # expected HALT3
    with MonkeyPatch.context() as m:
        mock_path = io.BytesIO(b"bad file contents")
        mock_path.exists = lambda: True
        test_component_monkey_patch = copy.deepcopy(test_component)
        test_component_monkey_patch.fastqGZ.path = mock_path
        flag = testCheck.validate(test_component_monkey_patch)
        assert flag.maxCode.name == "HALT3"
        assert (
            flag.message
            == "Corrupted Fastq.gz file suspected, last line number encountered: 0"
        )


def test_COMPONENT_RAWREADS_0001_single(
    glds48_dataSystem_STAGE00, test_COMPONENT_RAWREADS_0001
):
    ds = glds48_dataSystem_STAGE00
    testCheck = test_COMPONENT_RAWREADS_0001

    test_component = list(ds.dataset.samples.values())[0].rawReads

    # expected GREEN
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        flag = testCheck.validate(test_component)
        assert flag.maxCode.name == "GREEN"

    # expected HALT1
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        m.setattr(test_component, "fastqGZ", None)
        m.setattr(test_component, "fastqcReportZIP", None)
        flag = testCheck.validate(test_component)
        assert flag.maxCode.name == "HALT1"
        assert flag.message == "Missing expected files: ['fastqGZ', 'fastqcReportZIP']"

    # expected HALT2
    # TODO: requires a true truncated test file or something similar
    # with MonkeyPatch.context() as m:
    #    m.setattr(gzip, "open", lambda: raise EOFError)
    #    flag = COMPONENT_RAWREADS_0001.validate(test_component)
    #    assert flag.maxCode.name == 'HALT2'

    # expected HALT3
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        mock_path = io.BytesIO(b"bad file contents")
        mock_path.exists = lambda: True
        test_component_monkey_patch = copy.deepcopy(test_component)
        test_component_monkey_patch.fastqGZ.path = mock_path
        flag = testCheck.validate(test_component_monkey_patch)
        assert flag.maxCode.name == "HALT3"
        assert (
            flag.message
            == "Corrupted Fastq.gz file suspected, last line number encountered: 0"
        )


def test_COMPONENT_TRIMREADS_0001_paired(
    glds194_dataSystem_STAGE01, test_COMPONENT_TRIMREADS_0001
):
    ds = glds194_dataSystem_STAGE01
    testCheck = test_COMPONENT_TRIMREADS_0001

    test_component = list(ds.dataset.samples.values())[0].trimForwardReads

    # expected GREEN
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        flag = testCheck.validate(test_component)
        assert flag.maxCode.name == "GREEN"

    # expected HALT1
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        m.setattr(test_component, "trimmingReportTXT", None)
        flag = testCheck.validate(test_component)
        assert flag.maxCode.name == "HALT1"
        assert flag.message == "Missing expected files: ['trimmingReportTXT']"

    # expected HALT2
    # TODO: requires a true truncated test file or something similar
    # with MonkeyPatch.context() as m:
    #    m.setattr(gzip, "open", lambda: raise EOFError)
    #    flag = COMPONENT_RAWREADS_0001.validate(test_component)
    #    assert flag.maxCode.name == 'HALT2'

    # expected HALT3
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        mock_path = io.BytesIO(b"bad file contents")
        mock_path.exists = lambda: True
        test_component_monkey_patch = copy.deepcopy(test_component)
        test_component_monkey_patch.fastqGZ.path = mock_path
        flag = testCheck.validate(test_component_monkey_patch)
        assert flag.maxCode.name == "HALT3"
        assert (
            flag.message
            == "Corrupted Fastq.gz file suspected, last line number encountered: 0"
        )


def test_COMPONENT_TRIMREADS_0001_single(
    glds48_dataSystem_STAGE01, test_COMPONENT_TRIMREADS_0001
):
    ds = glds48_dataSystem_STAGE01
    testCheck = test_COMPONENT_TRIMREADS_0001

    test_component = list(ds.dataset.samples.values())[0].trimReads

    # expected GREEN
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        flag = testCheck.validate(test_component)
        assert flag.maxCode.name == "GREEN"

    # expected HALT1
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        m.setattr(test_component, "trimmingReportTXT", None)
        flag = testCheck.validate(test_component)
        assert flag.maxCode.name == "HALT1"
        assert flag.message == "Missing expected files: ['trimmingReportTXT']"

    # expected HALT2
    # TODO: requires a true truncated test file or something similar
    # with MonkeyPatch.context() as m:
    #    m.setattr(gzip, "open", lambda: raise EOFError)
    #    flag = COMPONENT_RAWREADS_0001.validate(test_component)
    #    assert flag.maxCode.name == 'HALT2'

    # expected HALT3
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        mock_path = io.BytesIO(b"bad file contents")
        mock_path.exists = lambda: True
        test_component_monkey_patch = copy.deepcopy(test_component)
        test_component_monkey_patch.fastqGZ.path = mock_path
        flag = testCheck.validate(test_component_monkey_patch)
        assert flag.maxCode.name == "HALT3"
        assert (
            flag.message
            == "Corrupted Fastq.gz file suspected, last line number encountered: 0"
        )


def test_COMPONENT_GENOMEALIGNMENTS_0001_paired(glds194_dataSystem_STAGE01):
    raise NotImplementedError


def test_COMPONENT_GENOMEALIGNMENTS_0001_single(glds48_dataSystem_STAGE01):
    raise NotImplementedError


def test_COMPONENT_RSEQCANALYSIS_0001_paired(glds194_dataSystem_STAGE01):
    raise NotImplementedError


def test_COMPONENT_RSEQCANALYSIS_0001_single(glds48_dataSystem_STAGE01):
    raise NotImplementedError


def test_COMPONENT_GENECOUNTS_0001_paired(glds194_dataSystem_STAGE01):
    raise NotImplementedError


def test_COMPONENT_GENECOUNTS_0001_single(glds48_dataSystem_STAGE01):
    raise NotImplementedError


def test_COMPONENT_DIFFERENTIALGENEEXPRESSION_0001_paired(glds194_dataSystem_STAGE01):
    raise NotImplementedError


def test_COMPONENT_DIFFERENTIALGENEEXPRESSION_0001_single(glds48_dataSystem_STAGE01):
    raise NotImplementedError

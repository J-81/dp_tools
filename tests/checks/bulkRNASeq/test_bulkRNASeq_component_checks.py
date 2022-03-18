import copy
import gzip
from io import StringIO
import io
import os
from pathlib import Path
from unittest.mock import MagicMock
import pickle
from dp_tools.bulkRNASeq.checks import COMPONENT_RAWREADS_0001, COMPONENT_TRIMREADS_0001
from dp_tools.bulkRNASeq.loaders import load_BulkRNASeq_STAGE_00
from dp_tools.core import check_model
from dp_tools.core.entity_model import DataFile

from pytest import MonkeyPatch


def test_COMPONENT_RAWREADS_0001_paired(glds194_dataSystem_STAGE00):
    ds = glds194_dataSystem_STAGE00

    test_component = list(ds.dataset.samples.values())[0].rawForwardReads

    # expected GREEN
    flag = COMPONENT_RAWREADS_0001.validate(test_component)
    assert flag.code.name == "GREEN"

    # expected HALT1
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        m.setattr(test_component, "fastqGZ", None)
        m.setattr(test_component, "fastqcReportZIP", None)
        flag = COMPONENT_RAWREADS_0001.validate(test_component)
        assert flag.code.name == "HALT1"
        assert flag.message == "Missing expected files: ['fastqGZ', 'fastqcReportZIP']"

    # expected HALT2
    # TODO: requires a true truncated test file or something similar
    # with MonkeyPatch.context() as m:
    #    m.setattr(gzip, "open", lambda: raise EOFError)
    #    flag = COMPONENT_RAWREADS_0001.validate(test_component)
    #    assert flag.code.name == 'HALT2'

    # expected HALT3
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        mock_path = io.BytesIO(b"bad file contents")
        mock_path.exists = lambda: True
        test_component_monkey_patch = copy.deepcopy(test_component)
        test_component_monkey_patch.fastqGZ.path = mock_path
        flag = COMPONENT_RAWREADS_0001.validate(test_component_monkey_patch)
        assert flag.code.name == "HALT3"
        assert (
            flag.message
            == "Corrupted Fastq.gz file suspected, last line number encountered: 0"
        )

def test_COMPONENT_RAWREADS_0001_single(glds48_dataSystem_STAGE00):
    ds = glds48_dataSystem_STAGE00

    test_component = list(ds.dataset.samples.values())[0].rawReads

    # expected GREEN
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        flag = COMPONENT_RAWREADS_0001.validate(test_component)
        assert flag.code.name == "GREEN"

    # expected HALT1
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        m.setattr(test_component, "fastqGZ", None)
        m.setattr(test_component, "fastqcReportZIP", None)
        flag = COMPONENT_RAWREADS_0001.validate(test_component)
        assert flag.code.name == "HALT1"
        assert flag.message == "Missing expected files: ['fastqGZ', 'fastqcReportZIP']"

    # expected HALT2
    # TODO: requires a true truncated test file or something similar
    # with MonkeyPatch.context() as m:
    #    m.setattr(gzip, "open", lambda: raise EOFError)
    #    flag = COMPONENT_RAWREADS_0001.validate(test_component)
    #    assert flag.code.name == 'HALT2'

    # expected HALT3
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        mock_path = io.BytesIO(b"bad file contents")
        mock_path.exists = lambda: True
        test_component_monkey_patch = copy.deepcopy(test_component)
        test_component_monkey_patch.fastqGZ.path = mock_path
        flag = COMPONENT_RAWREADS_0001.validate(test_component_monkey_patch)
        assert flag.code.name == "HALT3"
        assert (
            flag.message
            == "Corrupted Fastq.gz file suspected, last line number encountered: 0"
        )

def test_COMPONENT_TRIMREADS_0001_paired(glds194_dataSystem_STAGE01):
    ds = glds194_dataSystem_STAGE01

    test_component = list(ds.dataset.samples.values())[0].trimForwardReads

    # expected GREEN
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        flag = COMPONENT_TRIMREADS_0001.validate(test_component)
        assert flag.code.name == "GREEN"

    # expected HALT1
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        m.setattr(test_component, "trimmingReportTXT", None)
        flag = COMPONENT_TRIMREADS_0001.validate(test_component)
        assert flag.code.name == "HALT1"
        assert flag.message == "Missing expected files: ['trimmingReportTXT']"

    # expected HALT2
    # TODO: requires a true truncated test file or something similar
    # with MonkeyPatch.context() as m:
    #    m.setattr(gzip, "open", lambda: raise EOFError)
    #    flag = COMPONENT_RAWREADS_0001.validate(test_component)
    #    assert flag.code.name == 'HALT2'

    # expected HALT3
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        mock_path = io.BytesIO(b"bad file contents")
        mock_path.exists = lambda: True
        test_component_monkey_patch = copy.deepcopy(test_component)
        test_component_monkey_patch.fastqGZ.path = mock_path
        flag = COMPONENT_TRIMREADS_0001.validate(test_component_monkey_patch)
        assert flag.code.name == "HALT3"
        assert (
            flag.message
            == "Corrupted Fastq.gz file suspected, last line number encountered: 0"
        )



def test_COMPONENT_TRIMREADS_0001_single(glds48_dataSystem_STAGE01):
    ds = glds48_dataSystem_STAGE01

    test_component = list(ds.dataset.samples.values())[0].trimReads

    # expected GREEN
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        flag = COMPONENT_TRIMREADS_0001.validate(test_component)
        assert flag.code.name == "GREEN"

    # expected HALT1
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        m.setattr(test_component, "trimmingReportTXT", None)
        flag = COMPONENT_TRIMREADS_0001.validate(test_component)
        assert flag.code.name == "HALT1"
        assert flag.message == "Missing expected files: ['trimmingReportTXT']"

    # expected HALT2
    # TODO: requires a true truncated test file or something similar
    # with MonkeyPatch.context() as m:
    #    m.setattr(gzip, "open", lambda: raise EOFError)
    #    flag = COMPONENT_RAWREADS_0001.validate(test_component)
    #    assert flag.code.name == 'HALT2'

    # expected HALT3
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        mock_path = io.BytesIO(b"bad file contents")
        mock_path.exists = lambda: True
        test_component_monkey_patch = copy.deepcopy(test_component)
        test_component_monkey_patch.fastqGZ.path = mock_path
        flag = COMPONENT_TRIMREADS_0001.validate(test_component_monkey_patch)
        assert flag.code.name == "HALT3"
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
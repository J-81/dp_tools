""" Tests for validation report results, relies on test for loaders passing """
from pathlib import Path
import os

from pytest import MonkeyPatch
from dp_tools.bulkRNASeq.entity import BulkRNASeqSample

from dp_tools.bulkRNASeq.loaders import load_from_bulk_rnaseq_raw_dir
from dp_tools.bulkRNASeq.vv_protocols import BulkRNASeq_VVProtocol_RawData
from dp_tools.core.check_model import Check, Flag

# set for testing
TEST_DIR = Path(os.environ["TEST_ASSETS_DIR"])


def test_bulk_rnaseq_raw_data_validation(caplog):
    """ This tests validation as it would be run on dataset after demultiplexing """

    target_data_dir = TEST_DIR / "GLDS-194"

    caplog.set_level(0)
    ds = load_from_bulk_rnaseq_raw_dir(target_data_dir)
    vv_protocol = BulkRNASeq_VVProtocol_RawData(dataset=ds.dataset)

    with caplog.at_level(0):
        vv_protocol.validate_all()

    assert isinstance(vv_protocol.flags["dataset"], dict)
    assert isinstance(vv_protocol.flags["sample"], dict)
    assert isinstance(vv_protocol.flags["component"], dict)

    warn_msgs = [r.message for r in caplog.records if r.levelname == "WARNING"]
    # assert not warn_msgs  # no warnings expected

    # first, run a check with an unhandled developer exception
    with MonkeyPatch.context() as m:
        from dp_tools.bulkRNASeq.checks import SAMPLE_RAWREADS_0001

        m.setattr(SAMPLE_RAWREADS_0001, "validate_func", lambda: 1 / 0)
        vv_protocol.validate_all()
        sample, testcase_flags = list(vv_protocol.flags["sample"].items())[0]
        testcase_flag = testcase_flags[0]
        assert testcase_flag.check.id == "SAMPLE_RAWREADS_0001"
        assert testcase_flag.code.value == 91
        assert testcase_flag.message.startswith("An unhandled exception occured in ")
        assert isinstance(sample, BulkRNASeqSample)

    # second, run with only one check
    with MonkeyPatch.context() as m:
        vv_protocol.validate_all()
        sample, testcase_flags = list(vv_protocol.flags["sample"].items())[0]
        testcase_flag = testcase_flags[0]
        assert testcase_flag.check.id == "SAMPLE_RAWREADS_0001"
        assert testcase_flag.code.value == 20
        assert testcase_flag.message == "All expected raw read files present"
        assert isinstance(sample, BulkRNASeqSample)

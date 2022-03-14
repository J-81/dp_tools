""" Tests for specific validation checks """


from unittest.mock import MagicMock
from pytest import MonkeyPatch

from dp_tools.bulkRNASeq.checks import RAWREADS_0001
from dp_tools.bulkRNASeq.entity import BulkRNASeqSample
from dp_tools.components.components import ReadsComponent
from dp_tools.core import check_model


def test_RAWREADS_0001():
    # check attribute tests
    assert RAWREADS_0001.id == "RAWREADS_0001"

    # for paired end fails
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        mockSample = MagicMock(spec=BulkRNASeqSample)
        mockSample.dataset.metadata.paired_end = True
        flag = RAWREADS_0001.validate(sample=mockSample)
        assert flag.code == check_model.FlagCode.HALT1
        assert (
            flag.message
            == "Missing expected components: ['rawForwardReads', 'rawReverseReads']"
        )

    # for single end fails
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        mockSample = MagicMock(spec=BulkRNASeqSample)
        mockSample.dataset.metadata.paired_end = False
        flag = RAWREADS_0001.validate(sample=mockSample)
        assert flag.code == check_model.FlagCode.HALT1
        assert flag.message == "Missing expected components: ['rawReads']"

    # for paired end passes
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        mockSample = MagicMock(spec=BulkRNASeqSample)
        mockSample.rawForwardReads = MagicMock(spec=ReadsComponent)
        mockSample.rawReverseReads = MagicMock(spec=ReadsComponent)
        mockSample.dataset.metadata.paired_end = True
        flag = RAWREADS_0001.validate(sample=mockSample)
        assert flag.code == check_model.FlagCode.GREEN
        assert flag.message == "All expected raw read files present"

    # for single end passes
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        mockSample = MagicMock(spec=BulkRNASeqSample)
        mockSample.rawReads = MagicMock(spec=ReadsComponent)
        mockSample.dataset.metadata.paired_end = False
        flag = RAWREADS_0001.validate(sample=mockSample)
        assert flag.code == check_model.FlagCode.GREEN
        assert flag.message == "All expected raw read files present"

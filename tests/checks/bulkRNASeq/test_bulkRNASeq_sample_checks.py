""" Tests for specific validation checks """


from unittest.mock import MagicMock
from pytest import MonkeyPatch

from dp_tools.bulkRNASeq.checks import SAMPLE_RAWREADS_0001
from dp_tools.bulkRNASeq.entity import BulkRNASeqSample
from dp_tools.components.components import RawReadsComponent
from dp_tools.core import check_model


def test_SAMPLE_RAWREADS_0001():
    # check attribute tests
    assert SAMPLE_RAWREADS_0001.id == "SAMPLE_RAWREADS_0001"

    # for paired end fails
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        mockSample = MagicMock(spec=BulkRNASeqSample)
        mockForward = MagicMock(spec=RawReadsComponent)
        mockForward.count = 100
        mockReverse = MagicMock(spec=RawReadsComponent)
        mockReverse.count = 100
        mockSample.dataset.metadata.paired_end = True
        flag = SAMPLE_RAWREADS_0001.validate(sample=mockSample)
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
        flag = SAMPLE_RAWREADS_0001.validate(sample=mockSample)
        assert flag.code == check_model.FlagCode.HALT1
        assert flag.message == "Missing expected components: ['rawReads']"

    # for paired end passes
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        mockForward = MagicMock(spec=RawReadsComponent)
        mockForward.count = 100
        mockReverse = MagicMock(spec=RawReadsComponent)
        mockReverse.count = 100
        mockSample.rawForwardReads = mockForward
        mockSample.rawReverseReads = mockReverse
        mockSample.dataset.metadata.paired_end = True
        flag = SAMPLE_RAWREADS_0001.validate(sample=mockSample)
        assert flag.code == check_model.FlagCode.GREEN
        assert flag.message == "All expected raw read files present"

    # for paired end passes
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        mockForward = MagicMock(spec=RawReadsComponent)
        mockForward.count = 500
        mockReverse = MagicMock(spec=RawReadsComponent)
        mockReverse.count = 100
        mockSample.rawForwardReads = mockForward
        mockSample.rawReverseReads = mockReverse
        mockSample.dataset.metadata.paired_end = True
        flag = SAMPLE_RAWREADS_0001.validate(sample=mockSample)
        assert flag.code == check_model.FlagCode.HALT2
        assert flag.message == "Forward and reverse reads counts differ. Forward: (500) Reverse: (100)"

    # for single end passes
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        mockSample = MagicMock(spec=BulkRNASeqSample)
        mockSample.rawReads = MagicMock(spec=RawReadsComponent)
        mockSample.dataset.metadata.paired_end = False
        flag = SAMPLE_RAWREADS_0001.validate(sample=mockSample)
        assert flag.code == check_model.FlagCode.GREEN
        assert flag.message == "All expected raw read files present"

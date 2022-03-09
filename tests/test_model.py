from pathlib import Path
from dp_tools import model
import pytest
import logging
from unittest import mock

import os

# set for testing
TEST_DIR = Path(os.environ["TEST_ASSETS_DIR"])


def test_abc_inits_fails(caplog):
    caplog.set_level(logging.DEBUG)

    with pytest.raises(TypeError):
        model.BaseSample()


def test_empty_GLDSDataSystem():
    dSys = model.GLDSDataSystem(model.BaseDataSystem(name="GLDS-JJJ"))

    # empty sets expected
    assert dSys.all_datasets == set()
    assert dSys.all_samples == set()
    assert dSys.all_components == set()

    dSys.validate()


def test_datafile(caplog):
    target_data_file = (
        TEST_DIR
        / "GLDS-194"
        / "Metadata"
        / "AST_autogen_template_RNASeq_RCP_GLDS-194_RNASeq_runsheet.csv"
    )

    datf = model.DataFile(path=target_data_file)
    assert datf.md5sum == "ef85ec532b3e74b720131b3f5430443d"
    assert datf.path == target_data_file


def test_datafile_with_dummy_md5sum(caplog):
    target_data_file = (
        TEST_DIR
        / "GLDS-194"
        / "00-RawData"
        / "Fastq"
        / "Mmus_BAL-TAL_LRTN_BSL_Rep1_B7_R1_raw.fastq.gz"
    )

    datf = model.DataFile(path=target_data_file, dummy_md5sum=True)
    assert datf.md5sum == "DUMMY:8f3a03aed25aa00f091f3d4fd7fee3c0"
    assert datf.path == target_data_file


def test_bulk_rna_seq_dataset(caplog):
    dataset = model.BulkRNASeqDataset(base=model.BaseDataset(name="TestDSet"))

    bad_sample = "sample1"
    good_sample = mock.MagicMock(spec=model.BulkRNASeqSample)
    good_sample.name = "Sample2"

    # invalid sample: string
    with pytest.raises(TypeError):
        dataset.attach_sample(bad_sample)

    dataset.attach_sample(good_sample)

    # make sure a warning is logged if overwriting a sample
    with caplog.at_level(0):
        dataset.attach_sample(good_sample)
        assert caplog.records[-1].message == "Overwriting pre-existing sample: Sample2"

    dataset.validate()


def test_bulk_rna_seq_sample(caplog):
    caplog.set_level(0)
    sample = model.BulkRNASeqSample(base=model.BaseSample(name="test_sample_1"))

    assert sample.name == "test_sample_1"
    sample.validate()

    print("Before attach")
    sample.list_components()

    # looks like a component; however, doesn't have the base
    with pytest.raises(TypeError):
        mock_reads = mock.MagicMock(spec=model.ReadsComponent)
        sample.attach_component(mock_reads, attr="rawForwardReads")

    # looks like a component and has the proper base (at least a mock of one)
    mock_reads = mock.MagicMock(spec=model.ReadsComponent)
    mock_reads.base = mock.MagicMock(spec=model.BaseComponent)
    sample.attach_component(mock_reads, attr="rawForwardReads") 

    print("After attach")
    sample.list_components()
    1 / 0


def test_bulk_rna_seq_integration():
    dataset = model.BulkRNASeqDataset()
    sample = model.BulkRNASeqSample(base=model.BaseSample(name="test_sample_1"))

    dataset.attach_sample(sample)

    assert dataset.samples["test_sample_1"] == sample

    assert sample.dataset == dataset


def test_ingest_truncated_dataset():
    target_data_root = TEST_DIR / "GLDS-194"


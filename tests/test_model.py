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

    with pytest.raises(TypeError):
        model.BaseDataset()

    with pytest.raises(TypeError):
        model.BaseDataSystem()

    model.GLDSDataSystem(name="GLDS-JJJ")


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


def test_bulk_rna_seq_dataset():
    dataset = model.BulkRNASeqDataset()

    bad_sample = "sample1"
    good_sample = mock.MagicMock(spec=model.BulkRNASeqSample)
    good_sample.name = "Sample2"

    # invalid sample: string
    with pytest.raises(TypeError):
        dataset.attach_sample(bad_sample)

    dataset.attach_sample(good_sample)
    dataset.attach_sample(good_sample)

    dataset.validate()

    assert len(dataset.samples) == 1


def test_bulk_rna_seq_sample():
    sample = model.BulkRNASeqSample(base=model.BaseSample(name="test_sample_1"))

    assert sample.name == "test_sample_1"
    sample.validate()


def test_bulk_rna_seq_integration():
    dataset = model.BulkRNASeqDataset()
    sample = model.BulkRNASeqSample(base=model.BaseSample(name="test_sample_1"))

    dataset.attach_sample(sample)

    assert dataset.samples["test_sample_1"] == sample

    assert sample.dataset == dataset


def test_ingest_truncated_dataset():
    target_data_root = TEST_DIR / "GLDS-194"


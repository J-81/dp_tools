from dp_tools.bulkRNASeq.entity import BulkRNASeqDataset, BulkRNASeqSample
from dp_tools.components.components import RawReadsComponent
from dp_tools.core import entity_model
import pytest
from unittest import mock


def test_abc_init_fails():
    with pytest.raises(TypeError):
        entity_model.BaseSample()


def test_empty_GLDSDataSystem():
    dSys = entity_model.GLDSDataSystem(entity_model.BaseDataSystem(name="GLDS-JJJ"))

    # empty sets expected
    assert dSys.all_datasets == set()
    assert dSys.all_samples == set()
    assert dSys.all_components == set()

    # empty dict expected
    assert dSys.datasets == dict()


def test_datafile(glds194_runsheetPath):
    datf = entity_model.DataFile(path=glds194_runsheetPath)
    assert datf.compute_md5sum() == "02289bc0a96eb623c7850c07bcf7e196"
    assert datf.path == glds194_runsheetPath


def test_datafile_with_dummy_md5sum(glds194_runsheetPath):
    datf = entity_model.DataFile(path=glds194_runsheetPath, dummy_md5sum=True)
    assert datf.compute_md5sum() == "DUMMY:a5bfa1c912584b8a58e0fb9a56178778"
    assert datf.path == glds194_runsheetPath


def test_bulkRNASeq_dataset(caplog):
    dataset = BulkRNASeqDataset(base=entity_model.BaseDataset(name="TestDSet"))

    bad_sample = "sample1"
    good_sample = BulkRNASeqSample(base=mock.MagicMock(spec=entity_model.BaseSample))
    good_sample.base.name = "Sample2"

    # invalid sample: string
    with pytest.raises(TypeError):
        dataset.attach_sample(bad_sample)

    dataset.attach_sample(good_sample)

    # make sure a warning is logged if overwriting a sample
    with caplog.at_level(0):
        dataset.attach_sample(good_sample)
        assert caplog.records[-1].message == "Overwriting pre-existing sample: Sample2"


def test_bulkRNASeq_sample(caplog):
    caplog.set_level(0)
    sample = BulkRNASeqSample(base=entity_model.BaseSample(name="test_sample_1"))

    assert sample.name == "test_sample_1"
    # sample.validate()

    print("Before attach")
    assert isinstance(sample.rawForwardReads, entity_model.EmptyComponent)
    assert isinstance(sample.rawReverseReads, entity_model.EmptyComponent)

    # looks like a component; however, doesn't have the base
    with pytest.raises(TypeError):
        mock_reads = mock.MagicMock(spec=RawReadsComponent)
        sample.attach_component(mock_reads, attr="rawForwardReads")

    # looks like a component and has the proper base (at least a mock of one)
    mock_reads = mock.MagicMock(spec=RawReadsComponent)
    mock_reads.base = mock.MagicMock(spec=entity_model.BaseComponent)
    sample.attach_component(mock_reads, attr="rawForwardReads")

    print("After attach")
    assert isinstance(
        sample.rawForwardReads, RawReadsComponent
    )  # component should now be added
    assert isinstance(sample.rawReverseReads, entity_model.EmptyComponent)


def test_bulkRNASeq_full_noComponents(caplog):
    dataSystem = entity_model.GLDSDataSystem(
        base=entity_model.BaseDataSystem(name="Test_dataSystem_1")
    )
    dataset = BulkRNASeqDataset(base=entity_model.BaseDataset(name="Test_dataset_1"))
    sample = BulkRNASeqSample(base=entity_model.BaseSample(name="Test_sample_1"))

    dataset.attach_sample(sample)

    # before attachment this should return none
    assert dataSystem.dataset == None

    dataSystem.attach_dataset(dataset)

    # ensure warning logged during overwrite
    with caplog.at_level(0):
        dataSystem.attach_dataset(dataset)
        assert (
            caplog.records[-1].message
            == "Overwriting pre-existing dataset: Test_dataset_1__BulkRNASeq"
        )

    assert dataset.samples["Test_sample_1"] == sample

    assert sample.dataset == dataset

    # test all dataset accessors
    assert dataSystem.dataset == dataset
    assert dataSystem.datasets["Test_dataset_1__BulkRNASeq"] == dataset
    assert dataSystem.all_datasets == set([dataset])

    dataset2 = BulkRNASeqDataset(base=entity_model.BaseDataset(name="Test_dataset_2"))
    dataSystem.attach_dataset(dataset2)

    # test all dataset accessors
    # test all dataset accessors
    with pytest.raises(ValueError):
        dataSystem.dataset
    assert dataSystem.datasets["Test_dataset_1__BulkRNASeq"] == dataset
    assert dataSystem.all_datasets == set([dataset, dataset2])

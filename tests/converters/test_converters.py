from dp_tools.scripts.convert import isa_to_runsheet


def test_bulkRNASeq_STAGE00_validation_paired(glds194_dataSystem_STAGE00):
    """This tests validation as it would be run on dataset after demultiplexing"""
    ds = glds194_dataSystem_STAGE00
    df_runsheet = isa_to_runsheet(
        ds.name, ds.dataset.metadata.ISAarchive.path, config="0"
    )

    assert df_runsheet.shape == (13, 6)


def test_bulkRNASeq_STAGE00_validation_single(glds48_dataSystem_STAGE00):
    """This tests validation as it would be run on dataset after demultiplexing"""
    ds = glds48_dataSystem_STAGE00
    df_runsheet = isa_to_runsheet(
        ds.name, ds.dataset.metadata.ISAarchive.path, config="0"
    )

    assert df_runsheet.shape == (14, 6)

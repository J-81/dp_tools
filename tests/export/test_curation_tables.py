import os

from dp_tools.core.post_processing import update_curation_tables


def test_update_curation_tables_paired(glds194_dataSystem, tmpdir):
    os.chdir(tmpdir)
    df = update_curation_tables(
        glds194_dataSystem.dataset, config=("bulkRNASeq", "Latest")
    )

    assert df.shape == (13, 66)


def test_update_curation_tables_single(glds48_dataSystem, tmpdir):
    os.chdir(tmpdir)
    df = update_curation_tables(
        glds48_dataSystem.dataset, config=("bulkRNASeq", "Latest")
    )

    assert df.shape == (14, 56)

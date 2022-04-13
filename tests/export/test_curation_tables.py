from dp_tools.core.post_processing import update_curation_tables


def test_update_curation_tables(glds48_dataSystem_STAGE04):
    update_curation_tables(
        glds48_dataSystem_STAGE04.dataset, config=("bulkRNASeq", "0")
    )

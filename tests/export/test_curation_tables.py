from dp_tools.core.post_processing import export_curation_table


def test_export_curation_table(glds48_dataSystem_STAGE04):
    export_curation_table(glds48_dataSystem_STAGE04.dataset, config=("bulkRNASeq","0"))
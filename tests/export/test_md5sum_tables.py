from dp_tools.core.post_processing import (
    ALLOWED_MISSING_KEYS_FOR_NON_ERCC,
    ALLOWED_MISSING_KEYS_FOR_SINGLE_END,
    ALLOWED_MISSING_KEYS_FOR_PAIRED_END,
    generate_md5sum_table,
)


def test_generate_md5sum_tables_paired(glds194_dataSystem_STAGE04):
    missing_keys = ALLOWED_MISSING_KEYS_FOR_PAIRED_END
    df = generate_md5sum_table(
        glds194_dataSystem_STAGE04.dataset,
        config=("bulkRNASeq", "Latest"),
        allowed_unused_keys=missing_keys,
    )
    assert df.shape == (187, 3)


def test_generate_md5sum_tables_single(glds48_dataSystem_STAGE04):
    missing_keys = ALLOWED_MISSING_KEYS_FOR_NON_ERCC.union(
        ALLOWED_MISSING_KEYS_FOR_SINGLE_END
    )
    df = generate_md5sum_table(
        glds48_dataSystem_STAGE04.dataset,
        config=("bulkRNASeq", "Latest"),
        allowed_unused_keys=missing_keys,
    )

    assert df.shape == (153, 3)

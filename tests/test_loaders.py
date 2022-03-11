from pathlib import Path
import os

from dp_tools.bulkRNASeq.loaders import load_from_bulk_rnaseq_raw_dir

# set for testing
TEST_DIR = Path(os.environ["TEST_ASSETS_DIR"])


def test_from_bulk_rnaseq_raw_dir(caplog):
    target_data_dir = TEST_DIR / "GLDS-194"

    caplog.set_level(0)
    ds = load_from_bulk_rnaseq_raw_dir(
        target_data_dir
    )

    # pull dataset
    dataset = ds.datasets["GLDS-194__BulkRNASeq"]

    assert list(dataset.samples.keys()) == [
        "Mmus_BAL-TAL_LRTN_BSL_Rep1_B7",
        "Mmus_BAL-TAL_RRTN_BSL_Rep2_B8",
        "Mmus_BAL-TAL_RRTN_BSL_Rep3_B9",
        "Mmus_BAL-TAL_RRTN_BSL_Rep4_B10",
        "Mmus_BAL-TAL_LRTN_GC_Rep1_G6",
        "Mmus_BAL-TAL_LRTN_GC_Rep2_G8",
        "Mmus_BAL-TAL_LRTN_GC_Rep3_G9",
        "Mmus_BAL-TAL_RRTN_GC_Rep4_G10",
        "Mmus_BAL-TAL_LRTN_FLT_Rep1_F6",
        "Mmus_BAL-TAL_LRTN_FLT_Rep2_F7",
        "Mmus_BAL-TAL_LRTN_FLT_Rep3_F8",
        "Mmus_BAL-TAL_LRTN_FLT_Rep4_F9",
        "Mmus_BAL-TAL_LRTN_FLT_Rep5_F10",
    ]

    # check expected loaded components
    assert len(list(ds.dataset.samples["Mmus_BAL-TAL_LRTN_BSL_Rep1_B7"].rawForwardReads.multiQCDir.path.iterdir())) == 2
    assert len(list(ds.dataset.samples["Mmus_BAL-TAL_LRTN_BSL_Rep1_B7"].rawReverseReads.multiQCDir.path.iterdir())) == 2


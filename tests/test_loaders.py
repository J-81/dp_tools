from pathlib import Path
import os

from dp_tools.bulkRNASeq.loaders import load_from_bulk_rnaseq_raw_dir
import pytest

# set for testing
TEST_DIR = Path(os.environ["TEST_ASSETS_DIR"])

def test_from_bulk_rnaseq_with_bad_rootdir():
    target_data_dir = TEST_DIR / "GLDS-48_BUTWITHTYPOS"

    with pytest.raises(FileNotFoundError):
        ds = load_from_bulk_rnaseq_raw_dir(
            target_data_dir,
            dataSystem_name = "GLDS-48"
        )
    

def test_from_bulk_rnaseq_raw_dir_single(caplog):
    """ Tests loader for state after demultiplexing for single end study """
    target_data_dir = TEST_DIR / "GLDS-48_PostDemultiplex"

    caplog.set_level(0)
    ds = load_from_bulk_rnaseq_raw_dir(
        target_data_dir,
        dataSystem_name = "GLDS-48"
    )

    # pull dataset
    dataset = ds.datasets["GLDS-48__BulkRNASeq"]

    assert list(dataset.samples.keys()) == [
        'Mmus_C57-6J_LVR_GC_I_Rep1_M31', 'Mmus_C57-6J_LVR_GC_I_Rep2_M32', 'Mmus_C57-6J_LVR_FLT_I_Rep1_M21', 'Mmus_C57-6J_LVR_FLT_I_Rep2_M22', 'Mmus_C57-6J_LVR_GC_C_Rep1_M36', 'Mmus_C57-6J_LVR_GC_C_Rep2_M37', 'Mmus_C57-6J_LVR_GC_C_Rep3_M38', 'Mmus_C57-6J_LVR_GC_C_Rep4_M39', 'Mmus_C57-6J_LVR_GC_C_Rep5_M40', 'Mmus_C57-6J_LVR_FLT_C_Rep1_M25', 'Mmus_C57-6J_LVR_FLT_C_Rep2_M26', 'Mmus_C57-6J_LVR_FLT_C_Rep3_M27', 'Mmus_C57-6J_LVR_FLT_C_Rep4_M28', 'Mmus_C57-6J_LVR_FLT_C_Rep5_M30'
    ]

    # check expected loaded components [raw directory]
    assert len(list(ds.dataset.samples['Mmus_C57-6J_LVR_GC_I_Rep1_M31'].rawReads.multiQCDir.path.iterdir())) == 2


def test_from_bulk_rnaseq_raw_dir_paired(caplog):
    """ Tests loader for state after demultiplexing for single end study """
    target_data_dir = TEST_DIR / "GLDS-207_PostDemultiplex"

    caplog.set_level(0)
    ds = load_from_bulk_rnaseq_raw_dir(
        target_data_dir,
        dataSystem_name = "GLDS-207"
    )

    # pull dataset
    dataset = ds.datasets["GLDS-207__BulkRNASeq"]

    assert list(dataset.samples.keys()) == [
        'KCNQ370_Female_Ground_Control-01', 'KCNQ370_Female_Ground_Control-02', 'KCNQ370_Female_Ground_Control-03', 'KCNQ370_Female_Ground_Control-04', 'KCNQ97_Female_Ground_Control-05', 'KCNQ97_Female_Ground_Control-06', 'KCNQ97_Female_Ground_Control-07', 'KCNQ97_Female_Ground_Control-08', 'Sei_ts1-Female_Ground_Control-09', 'Sei_ts1-Female_Ground_Control-10', 'Sei_ts1-Female_Ground_Control-11', 'Sei_ts1-Female_Ground_Control-12', 'CS-Female_Ground_Control-13', 'CS-Female_Ground_Control-14', 'CS-Female_Ground_Control-15', 'CS-Female_Ground_Control-16', 'KCNQ370_Female_Space_Flown-17', 'KCNQ370_Female_Space_Flown-18', 'KCNQ370_Female_Space_Flown-19', 'KCNQ370_Female_Space_Flown-20', 'KCNQ97_Female_Space_Flown-21', 'KCNQ97_Female_Space_Flown-22', 'KCNQ97_Female_Space_Flown-23', 'KCNQ97_Female_Space_Flown-24', 'Sei_ts1-Female_Space_Flown-25', 'Sei_ts1-Female_Space_Flown-26', 'Sei_ts1-Female_Space_Flown-27', 'Sei_ts1-Female_Space_Flown-28', 'CS-Female_Space_Flown-29', 'CS-Female_Space_Flown-30', 'CS-Female_Space_Flown-31', 'CS-Female_Space_Flown-32'
    ]

    # check expected loaded components [raw directory]
    for sample_name in ds.dataset.samples:
        assert len(list(ds.dataset.samples[sample_name].rawForwardReads.multiQCDir.path.iterdir())) == 2
        assert len(list(ds.dataset.samples[sample_name].rawReverseReads.multiQCDir.path.iterdir())) == 2
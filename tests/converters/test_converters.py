import hashlib

import pandas as pd

from dp_tools.scripts.convert import isa_to_runsheet

# Updated to microarray development version api
def test_paired_isa_to_runsheet(glds194_dataSystem_STAGE00):
    """This tests validation as it would be run on dataset after demultiplexing"""
    ds = glds194_dataSystem_STAGE00
    df_runsheet = isa_to_runsheet(
        ds.name, ds.dataset.metadata.ISAarchive.path, config=("bulkRNASeq", "0")
    )

    assert df_runsheet.shape == (13, 7)
    assert (
        hashlib.sha1(pd.util.hash_pandas_object(df_runsheet).values).hexdigest()
        == "99be45a88523dd07ac17e18ab13838ef7a24e883"
    ), "Hash did not match, the means the contents changed. Manually validation and reset of test hash is in order"


def test_single_isa_to_runsheet(glds48_dataSystem_STAGE00):
    """This tests validation as it would be run on dataset after demultiplexing"""
    ds = glds48_dataSystem_STAGE00
    df_runsheet = isa_to_runsheet(
        ds.name, ds.dataset.metadata.ISAarchive.path, config=("bulkRNASeq", "0")
    )

    assert df_runsheet.shape == (14, 7)
    assert (
        hashlib.sha1(pd.util.hash_pandas_object(df_runsheet).values).hexdigest()
        == "0477c47d0b081e5a944d6f343d1737bdfa8dd274"
    ), "Hash did not match, the means the contents changed. Manually validation and reset of test hash is in order"


def test_microarray_glds205_isa_to_runsheet(glds205_isazip_path):
    """This tests validation as it would be run on dataset after demultiplexing"""
    df_runsheet = isa_to_runsheet(
        "GLDS-205", glds205_isazip_path, config=("microarray", "0")
    )

    assert df_runsheet.shape == (16, 13)  # 2 factor values
    assert (
        hashlib.sha1(pd.util.hash_pandas_object(df_runsheet).values).hexdigest()
        == "f8bab4078c8b4715708f6fccc1c314d853f4fb9c"
    ), "Hash did not match, the means the contents changed. Manually validation and reset of test hash is in order"


def test_microarray_glds123_isa_to_runsheet(glds123_isazip_path):
    """This tests validation as it would be run on dataset after demultiplexing"""
    df_runsheet = isa_to_runsheet(
        "GLDS-123", glds123_isazip_path, config=("microarray", "0")
    )

    assert df_runsheet.shape == (16, 12)  # 1 factor values
    assert (
        hashlib.sha1(pd.util.hash_pandas_object(df_runsheet).values).hexdigest()
        == "2cd5f640d5b0ced45be1f845e0aa9313bd2e0803"
    ), "Hash did not match, the means the contents changed. Manually validation and reset of test hash is in order"

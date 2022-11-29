import hashlib
import os

import pandas as pd

from dp_tools.scripts.convert import isa_to_runsheet

# Updated to microarray development version api
def test_paired_isa_to_runsheet(glds194_test_dir, tmpdir):
    """This tests validation as it would be run on dataset after demultiplexing"""
    os.chdir(tmpdir)
    isaPath = glds194_test_dir / "Metadata" / "GLDS-194_metadata_GLDS-194-ISA.zip"
    df_runsheet = isa_to_runsheet("GLDS-194", isaPath, config=("bulkRNASeq", "0"))

    assert df_runsheet.shape == (13, 8)
    assert (
        hashlib.sha1(pd.util.hash_pandas_object(df_runsheet).values).hexdigest()
        == "4b4368ebf60c9e02a11b792030102f8d67e43b67"
    ), "Hash did not match, the means the contents changed. Manually validation and reset of test hash is in order"


def test_single_isa_to_runsheet(glds48_test_dir, tmpdir):
    """This tests validation as it would be run on dataset after demultiplexing"""
    os.chdir(tmpdir)
    isaPath = glds48_test_dir / "Metadata" / "GLDS-48_metadata_RR1-NASA-ISA.zip"
    df_runsheet = isa_to_runsheet("GLDS-48", isaPath, config=("bulkRNASeq", "0"))

    assert df_runsheet.shape == (14, 7)
    assert (
        hashlib.sha1(pd.util.hash_pandas_object(df_runsheet).values).hexdigest()
        == "8322bdba7204e3f685b4af8affec5da9fc5bd526"
    ), "Hash did not match, the means the contents changed. Manually validation and reset of test hash is in order"


def test_methylSeq_glds397_isa_to_runsheet(glds397_isazip_path):
    """This tests isa_to_runsheet on a methylSeq assay, GLDS-397"""
    df_runsheet = isa_to_runsheet(
        "GLDS-397", glds397_isazip_path, config = ("methylSeq", "1")
    )

    assert df_runsheet.shape == (16, 6)  # 1 factor value
    assert (
        hashlib.sha1(pd.util.hash_pandas_object(df_runsheet).values).hexdigest()
        == "57cab80a811540e3f1174ac6b65c7ecc3e7f73f5"
    ), "Hash did not match, the means the contents changed. Manually validation and reset of test hash is in order"



# TODO: Enable once microarray test data is included in test repo
# def test_microarray_glds205_isa_to_runsheet(glds205_isazip_path):
#     """This tests validation as it would be run on dataset after demultiplexing"""
#     df_runsheet = isa_to_runsheet(
#         "GLDS-205", glds205_isazip_path, config=("microarray", "0")
#     )

#     assert df_runsheet.shape == (16, 13)  # 2 factor values
#     assert (
#         hashlib.sha1(pd.util.hash_pandas_object(df_runsheet).values).hexdigest()
#         == "c316e9e8b6013795395fab270a2bc96bd917d778"
#     ), "Hash did not match, the means the contents changed. Manually validation and reset of test hash is in order"


# def test_microarray_glds123_isa_to_runsheet(glds123_isazip_path):
#     """This tests validation as it would be run on dataset after demultiplexing"""
#     df_runsheet = isa_to_runsheet(
#         "GLDS-123", glds123_isazip_path, config=("microarray", "0")
#     )

#     assert df_runsheet.shape == (16, 12)  # 1 factor values
#     assert (
#         hashlib.sha1(pd.util.hash_pandas_object(df_runsheet).values).hexdigest()
#         == "0dce69a3de3783856f7f6acc7f1eb9e634d4e2c8"
#     ), "Hash did not match, the means the contents changed. Manually validation and reset of test hash is in order"

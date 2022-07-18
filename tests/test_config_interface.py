from dp_tools.config.interface import *

def test_load_config():
    res = load_config(("bulkRNASeq","0"))
    assert res["data assets"]
    res = load_config(("bulkRNASeq","1"))
    assert res["data assets"]
    res = load_config(("bulkRNASeq","Latest"))
    assert res["data assets"]

    #res = load_config(("microarray","0"))
    #assert res["data assets"]

def test_get_data_asset_template():
    ...

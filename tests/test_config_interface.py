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

def test_get_data_asset_keys():
    res = get_data_asset_keys(config = ("bulkRNASeq","Latest"))
    assert 'raw reads fastq GZ' in res

def test_get_data_asset_template():
    res = get_data_asset_template(key = "raw reads fastq GZ",config = ("bulkRNASeq","Latest"))
    assert res == "00-RawData/Fastq/{sample}_raw.fastq.gz"
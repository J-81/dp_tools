from dp_tools.glds_api.files import get_urls


def test_get_full_filelisting():
    accession = "GLDS-194"
    files = get_urls(accession)
    assert len(files) == 289

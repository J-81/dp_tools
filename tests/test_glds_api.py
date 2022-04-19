from dp_tools.glds_api.files import associate_url_for_data_file, get_urls
import pytest


def test_get_full_filelisting():
    accession = "GLDS-194"
    files = get_urls(accession)


def test_get_url_for_data_file(glds194_dataSystem_STAGE00):
    ds = glds194_dataSystem_STAGE00.dataset
    for sample in ds.samples.values():
        associate_url_for_data_file(sample, sample.rawForwardReads.fastqGZ)
        assert sample.rawForwardReads.fastqGZ.glds_url

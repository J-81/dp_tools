import pytest

from dp_tools.glds_api.commons import (
    get_table_of_files,
    retrieve_filenames_for_datatype,
    retrieve_datatypes_from_accession,
    retrieve_file_url,
)


def test_get_table_of_files():
    accession = "GLDS-194"
    df = get_table_of_files(accession)
    assert len(df) == 221

    accession = "GLDS-1"
    df = get_table_of_files(accession)
    assert len(df) == 144


def test_retrieve_filenames_for_datatype():
    accession = "GLDS-194"
    filenames = retrieve_filenames_for_datatype(accession, datatype="isa")
    assert len(filenames) == 1

    accession = "GLDS-194"
    filenames = retrieve_filenames_for_datatype(accession, datatype="raw reads")
    assert len(filenames) == 26

    accession = "GLDS-1"
    filenames = retrieve_filenames_for_datatype(accession, datatype="no_type")
    assert len(filenames) == 0


def test_retrieve_datatypes_from_accession():
    accession = "GLDS-194"
    datatypes = retrieve_datatypes_from_accession(accession)
    assert len(datatypes) == 14

    accession = "GLDS-1"
    datatypes = retrieve_datatypes_from_accession(accession)
    assert len(datatypes) == 8


def test_retrieve_file_url():
    accession = "GLDS-194"
    url = retrieve_file_url(accession, filename="GLDS-194_metadata_GLDS-194-ISA.zip")
    assert (
        url
        == "https://genelab-data.ndc.nasa.gov/geode-py/ws/studies/GLDS-194/download?file=GLDS-194_metadata_GLDS-194-ISA.zip"
    )

    accession = "GLDS-1"
    # Use a non-existent filename to trigger exception
    with pytest.raises(ValueError):
        url = retrieve_file_url(
            accession, filename="GLDS-194_metadata_GLDS-194-ISA.zip"
        )

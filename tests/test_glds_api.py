import pytest

from dp_tools.glds_api.commons import (
    get_table_of_files,
    retrieve_file_url,
    find_matching_filenames
)


def test_get_table_of_files():
    accession = "GLDS-194"
    df = get_table_of_files(accession)
    assert len(df) == 285

    accession = "GLDS-1"
    df = get_table_of_files(accession)
    assert len(df) == 72

def test_retrieve_file_url():
    accession = "GLDS-194"
    url = retrieve_file_url(accession, filename="GLDS-194_metadata_GLDS-194-ISA.zip")
    assert (
        url
        == "https://genelab-data.ndc.nasa.gov/geode-py/ws/studies/OSD-194/download?source=datamanager&file=GLDS-194_metadata_GLDS-194-ISA.zip"
    )

    accession = "GLDS-1"
    # Use a non-existent filename to trigger exception
    with pytest.raises(ValueError):
        url = retrieve_file_url(
            accession, filename="GLDS-194_metadata_GLDS-194-ISA.zip"
        )

def test_find_matching_filenames():
    accession = "GLDS-194"
    filenames = find_matching_filenames(accession, filename_pattern=".*-ISA.zip")
    assert (
        filenames
        == ['GLDS-194_metadata_GLDS-194-ISA.zip']
    )

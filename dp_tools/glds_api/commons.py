"""
Python functions the retrieve data from GeneLab. Uses the GeneLab openAPI (https://visualization.genelab.nasa.gov/GLOpenAPI) to retrieve most data

Direct file retrieval is supported using URLs present on the GLDS pages themselves.
"""

import functools
from urllib.request import urlopen
import logging

import pandas as pd

log = logging.getLogger(__name__)

GENELAB_DATASET_FILES = "https://visualization.genelab.nasa.gov/GLOpenAPI/samples/?id={accession}&file.datatype&format=csv"
""" Template URL to access table of filenames for a single GLDS accession ID """

FILE_RETRIEVAL_URL = "https://genelab-data.ndc.nasa.gov/geode-py/ws/studies/{accession}/download?file={filename}"
""" Template URL to download a certain file by filename

Note: This url is not part of the openAPI and will be replaced in the future.
"""


@functools.cache
def get_table_of_files(accession: str) -> pd.DataFrame:
    """Retrieve table of filenames associated with a GLDS accesion ID.

    Note: This function is cached to prevent extra api calls. This can desync from the GLDS repository in the rare case that the GLDS accession is updated in between related calls.

    :param accession: GLDS accession ID, e.g. 'GLDS-194'
    :type accession: str
    :return: A dataframe containing each filename including associated metadata like datatype
    :rtype: pd.DataFrame
    """
    # format url
    log.info(f"Retrieving table of files for {accession}")
    url = GENELAB_DATASET_FILES.format(accession=accession)

    # fetch data
    log.info(f"URL Source: {url}")
    with urlopen(url) as response:
        df = pd.read_csv(response, header=1)
    return df


def retrieve_file_url(accession: str, filename: str) -> str:
    """Retrieve file URL associated with a GLDS accesion ID

    :param accession: GLDS accession ID, e.g. 'GLDS-194'
    :type accession: str
    :param filename: Full filename, e.g. 'GLDS-194_metadata_GLDS-194-ISA.zip'
    :type filename: str
    :return: URL to fetch the most recent version of the file
    :rtype: str
    """
    # Check that the filenames exists
    df = get_table_of_files(accession)
    if filename not in list(df["filename"]):
        raise ValueError(
            f"Could not find filename: '{filename}'. Here as are found filenames for '{accession}': '{df['filename'].unique()}'"
        )

    return FILE_RETRIEVAL_URL.format(accession=accession, filename=filename)


def retrieve_filenames_for_datatype(accession: str, datatype: str) -> list[str]:
    """Retrieve filenames for a specific datatype and accession ID

    :param accession: GLDS accession ID, e.g. 'GLDS-194'
    :type accession: str
    :param datatype: datatype of the files requested, e.g. 'isa'
    :type datatype: str
    :return: URL to fetch the most recent version of the file
    :rtype: list[str]
    """
    df = get_table_of_files(accession)

    # filter to only requested datatype
    df = df.loc[df["datatype"] == datatype]

    # convert filename column into list of unique filenames
    filenames = list(df["filename"].unique())

    return filenames


def retrieve_datatypes_from_accession(accession: str) -> list[str]:
    """Retrieve all datatypes for an accession ID

    :param accession: GLDS accession ID, e.g. 'GLDS-194'
    :type accession: str
    :return: Data types associated with the accesssion ID
    :rtype: list[str]
    """
    df = get_table_of_files(accession)

    # convert filename column into list of unique filenames
    datatypes = list(df["datatype"].unique())

    return datatypes

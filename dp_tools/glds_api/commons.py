"""
Python functions the retrieve data from GeneLab. Uses the GeneLab public APIs (https://genelab.nasa.gov/genelabAPIs)
"""

import functools
from urllib.request import urlopen
import logging

import yaml
import pandas as pd

log = logging.getLogger(__name__)

GENELAB_DATASET_FILES = "https://genelab-data.ndc.nasa.gov/genelab/data/glds/files/{accession_number}"
""" Template URL to access json of files for a single GLDS accession ID """

FILE_RETRIEVAL_URL_PREFIX = "https://genelab-data.ndc.nasa.gov{suffix}"
""" Used to retrieve files using remote url suffixes listed in the 'Data Query' API """

@functools.cache
def get_table_of_files(accession: str) -> pd.DataFrame:
    """Retrieve table of filenames associated with a GLDS accesion ID.
    Converts 'GLDS' to 'OSD' as required by open sciences repository change

    Note: This function is cached to prevent extra api calls. This can desync from the GLDS repository in the rare case that the GLDS accession is updated in between related calls.

    :param accession: GLDS accession ID, e.g. 'GLDS-194'
    :type accession: str
    :return: A dataframe containing each filename including associated metadata like datatype
    :rtype: pd.DataFrame
    """
    # format url
    log.info(f"Retrieving table of files for {accession}")
    url = GENELAB_DATASET_FILES.format(accession_number=accession.split("-")[1])

    accession_osd = accession.replace("GLDS", "OSD") # this is the new study level accession ID

    # fetch data
    log.info(f"URL Source: {url}")
    with urlopen(url) as response:
        data = yaml.safe_load(response.read())
        df = pd.DataFrame(data['studies'][accession_osd]['study_files'])
    return df

def find_matching_filenames(accession: str, filename_pattern: str) -> list[str]:
    """Returns list of file names that match the provided regex pattern.

    :param accession: GLDS accession ID, e.g. 'GLDS-194'
    :type accession: str
    :param filename_pattern: Regex pattern to query against file names
    :type filename_pattern: str
    :return: List of file names that match the regex
    :rtype: list[str]
    """
    df = get_table_of_files(accession)
    return df.loc[df['file_name'].str.contains(filename_pattern), 'file_name'].to_list()

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
    if filename not in list(df["file_name"]):
        raise ValueError(
            f"Could not find filename: '{filename}'. Here as are found filenames for '{accession}': '{df['file_name'].unique()}'"
        )
    url = FILE_RETRIEVAL_URL_PREFIX.format(suffix=df.loc[df['file_name'] == filename, 'remote_url'].squeeze())
    print(url)
    return url

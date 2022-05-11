import datetime
import functools
from json import loads
from urllib.request import urlopen
import logging

log = logging.getLogger(__name__)
# Function to pull metadata zip from GeneLab
# Credit to Kirill Grigorev
GENELAB_ROOT = "https://genelab-data.ndc.nasa.gov"
GLDS_URL_PREFIX = GENELAB_ROOT + "/genelab/data/study/data/"
FILELISTINGS_URL_PREFIX = GENELAB_ROOT + "/genelab/data/study/filelistings/"
ISA_ZIP_REGEX = r".*_metadata_.*[_-]ISA\.zip$"


def read_json(url):
    with urlopen(url) as response:
        return loads(response.read().decode())


@functools.cache
def get_glds_filelisting_json(accession: str) -> tuple[dict, str]:
    """Return filelisting json accession number"""
    log.debug(f"Fetching filelisting JSON for '{accession}'")
    # extract file urls
    url1 = GLDS_URL_PREFIX + accession
    log.debug(f"Using {url1} to find study internal id for {accession}")
    glds_json = read_json(url1)
    try:
        _id = glds_json[0]["_id"]
    except (AssertionError, TypeError, KeyError, IndexError):
        raise ValueError("Malformed JSON?")
    file_list_url = FILELISTINGS_URL_PREFIX + _id
    file_listing_json = read_json(file_list_url)
    version = glds_json[0]["version"]
    log.info(
        f"Retrieved file listing json. URL: '{file_list_url}' Version: '{version}' Access Date: '{datetime.datetime.now().strftime('%Y-%m-%d %H:%M %Z')}'"
    )
    return file_listing_json

import argparse
from pathlib import Path
from re import search
from urllib.parse import quote
from dp_tools.glds_api.commons import (
    FILELISTINGS_URL_PREFIX,
    GENELAB_ROOT,
    GLDS_URL_PREFIX,
    ISA_ZIP_REGEX,
    read_json,
)

import requests
import logging

log = logging.getLogger(__name__)


def _parse_args():
    """Parse command line args."""
    parser = argparse.ArgumentParser(
        description=f"Script for downloading latest ISA from GLDS repository"
    )
    parser.add_argument(
        "--accession", metavar="GLDS-001", required=True, help="GLDS accession number"
    )
    parser.add_argument(
        "--alternate-url",
        action="store_true",
        default=False,
        help="Use alternate url, fetched by api script",
    )

    args = parser.parse_args()
    return args


def get_isa(accession: str) -> tuple[str, str, str, str]:
    """Returns isa filename as well as GeneLab URLS from the associated file listing

    :param accession: GLDS accession ID, e.g. GLDS-194
    :type accession: str
    :raises ValueError: "Malformed JSON?" if there is an issue getting the filelisting url 'id'
    :raises ValueError: "Unexpected: no ISAs found"
    :raises ValueError: "Unexpected: multiple files match the ISA regex"
    :return: A tuple as follows: [filename, version_id, url, alternative_url]
    :rtype: tuple[str, str, str, str]
    """
    glds_json = read_json(GLDS_URL_PREFIX + accession)
    try:
        _id = glds_json[0]["_id"]
        log.info(f"File Listing URL: {FILELISTINGS_URL_PREFIX + _id}")
    except (AssertionError, TypeError, KeyError, IndexError):
        raise ValueError("Malformed JSON?")
    isa_entries = [
        entry
        for entry in read_json(FILELISTINGS_URL_PREFIX + _id)
        if search(ISA_ZIP_REGEX, entry["file_name"])
    ]
    if len(isa_entries) == 0:
        raise ValueError("Unexpected: no ISAs found")
    elif len(isa_entries) > 1:
        raise ValueError("Unexpected: multiple files match the ISA regex")
    else:
        entry = isa_entries[0]
        version = entry["version"]
        url = GENELAB_ROOT + entry["remote_url"] + "?version={}".format(version)
        alt_url = (
            GENELAB_ROOT
            + "/genelab/static/media/dataset/"
            + quote(entry["file_name"])
            + "?version={}".format(version)
        )
        return entry["file_name"], version, url, alt_url


def download_isa(accession: str, alternate_url: bool = False) -> str:
    """Downloads the ISA archive for the given GLDS accession number.

    :param accession: GLDS accession number, e.g. GLDS-194
    :type accession: str
    :param alternate_url: Designates if the alternative url should be used, both fetch the same file, defaults to False
    :type alternate_url: bool, optional
    :return: The filename as downloaded to the local drive
    :rtype: str
    """

    log.info(f"Accessing GeneLab API for ISA file. Accession: {accession}")
    filename, _, url, alt_url = get_isa(accession)
    if not Path(filename).is_file():
        log.info(f"Successfully retrieved ISA file location from API.")
        use_url = url if not alternate_url else alt_url
        log.info(f"Downloading from {use_url}.")
        r = requests.get(use_url)
        # If the response was successful, no Exception will be raised
        r.raise_for_status()
        with open(filename, "wb") as f:
            f.write(r.content)
        log.info(f"Finished downloading ISA file: {filename}")
    else:
        log.info(f"Already downloaded {filename} to current folder")
    return filename


def main():
    args = _parse_args()
    download_isa(args.accession, args.alternate_url)


if __name__ == "__main__":
    isazip = main()

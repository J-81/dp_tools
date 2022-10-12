import argparse
from pathlib import Path
from dp_tools.glds_api.commons import retrieve_file_url, find_matching_filenames

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

    args = parser.parse_args()
    return args


def download_isa(accession: str) -> str:
    """Downloads the ISA archive for the given GLDS accession number.

    :param accession: GLDS accession number, e.g. GLDS-194
    :type accession: str
    :param alternate_url: Designates if the alternative url should be used, both fetch the same file, defaults to False
    :type alternate_url: bool, optional
    :return: The filename as downloaded to the local drive
    :rtype: str
    """

    log.info(f"Accessing GeneLab API for ISA file. Accession: {accession}")
    [filename] = find_matching_filenames(accession, filename_pattern = ".*-ISA.zip")
    url = retrieve_file_url(accession, filename = filename)
    if not Path(filename).is_file():
        log.info(f"Successfully retrieved ISA file location from API.")
        log.info(f"Downloading from {url}.")
        r = requests.get(url)
        # If the response was successful, no Exception will be raised
        r.raise_for_status()
        with open(filename, "wb") as f:
            f.write(r.content)
        log.info(f"Finished downloading ISA file: {filename}")
    else:
        log.info(f"Already downloaded {filename} to current folder")
    return filename


def main():
    logging.basicConfig(level=logging.INFO)
    args = _parse_args()
    download_isa(args.accession)


if __name__ == "__main__":
    isazip = main()

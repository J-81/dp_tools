import sys

import requests
import click

from dp_tools.glds_api.commons import find_matching_filenames, retrieve_file_url

@click.command()
@click.option("--osd-id", help='OSD Accession ID. e.g. "OSD-194"', required=True)
def main(osd_id):
    files = find_matching_filenames(accession=osd_id, filename_pattern=".*multiqc.*.zip")
    # Ensure we also download ISA archive
    files.extend(
        find_matching_filenames(accession=osd_id, filename_pattern=".*ISA.*.zip")
    )

    if not any(["read_dist" in f for f in files]):
        print(
            "Did not locate read distribution RSeQC report zip.  Inferring this isn't a dataset to be used."
        )
        sys.exit(0)


    def download_file(url, local_filename):
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(local_filename, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)


    for f in files:
        download_file(retrieve_file_url(osd_id, f), f)

if __name__ == '__main__':
    main()
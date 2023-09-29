import sys
from pathlib import Path
import zipfile

import requests
import click

from dp_tools.glds_api.commons import find_matching_filenames, retrieve_file_url

@click.command()
@click.option("--osd-id", help='OSD Accession ID. e.g. "OSD-194"', required=True)
@click.option("--output-dir", help='Output directory', required=False, default=".")
def main(osd_id, output_dir):


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
        print(f"Saving file: {local_filename} from {url}")
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(local_filename, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)

    def unzip_file(zip_file_name, extraction_location):
        print(f"Unzipping file: {zip_file_name}")
        # open the zip file in read mode
        with zipfile.ZipFile(zip_file_name, 'r') as zip_ref:
            # extract all the contents into the directory
            zip_ref.extractall(extraction_location)

    # Setup output dir
    output_dir = Path(output_dir)
    if not output_dir.exists():
        output_dir.mkdir()
    for f in files:
        file_location = output_dir / f
        download_file(retrieve_file_url(osd_id, f), file_location)
        unzip_file(file_location, output_dir)

if __name__ == '__main__':
    main()
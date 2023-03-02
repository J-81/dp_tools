import enum
from pathlib import Path
import sys
import requests

import click
from loguru import logger
import pandas as pd

from dp_tools.glds_api import commons, isa
from dp_tools.core.files import isa_archive

@click.group()
def osd():
    pass


@click.command()
@click.argument("osd-id")
@click.argument("file-pattern")
@click.option("--dry-run", default=False, is_flag=True, help="Output urls and filenames only.")
@click.option("-y","--y","non_interactive", default=False, is_flag=True, help="Download files without a prompt")
def download_files(osd_id, file_pattern, dry_run, non_interactive):
    
    # Determine table of filenames
    logger.info(f"Fetching file list for {osd_id}")

    filenames = commons.find_matching_filenames(accession = osd_id, filename_pattern=file_pattern)

    logger.info(f"Found {len(filenames)} associated with {osd_id} after filtering based on regex pattern: {file_pattern}")
    
    # For file, get url and echo
    if dry_run:
        logger.info(f"Output pairs of 'url' 'filename'")
        for filename in filenames:
            url = commons.retrieve_file_url(accession = osd_id, filename = filename)
            click.echo(f"{url} {filename}")
    else:
        if not non_interactive:
            click.prompt("Ready to download files... Press enter to continue.")
        for i, filename in enumerate(filenames):
            url = commons.retrieve_file_url(accession = osd_id, filename = filename)
            logger.info(f"Downloading file: {filename}. {i+1} of {len(filenames)}")
            logger.debug(f"Download url: {url}")
            r = requests.get(url)  
            with open(filename, 'wb') as f:
                f.write(r.content)
@click.command()
@click.argument("osd-id")
@click.option("-o", "--output", default="samples.txt", show_default=True, help="File to write sample names.")
def get_samples(osd_id, output):
    logger.info(f"Fetching information for {osd_id}")

    isa_path = isa.download_isa(accession = osd_id)

    isa_sample_and_assay_files = [f for f in isa_archive.fetch_isa_files(Path(isa_path)) if any([f.name.startswith("a_"), f.name.startswith("s_")])]
    isa_sample_and_assay_files.sort()

    logger.info(f"Found these ISA assay files: {[f.name for f in isa_sample_and_assay_files]}")

    isa_sample_and_assay_files = {i:f for i,f in enumerate(isa_sample_and_assay_files)}

    for i,f in isa_sample_and_assay_files.items():
        click.echo(f"{i}: {f.name}", err=True)
    selection = click.prompt("Select an table file by number", type=int, err=True)
    if selection not in isa_sample_and_assay_files:
        raise ValueError("Invalid choice!")
    else:
        target_file = isa_sample_and_assay_files[selection]
        logger.info(f"Selected {target_file}")
    
    samples = [s.strip() for s in pd.read_csv(target_file, sep="\t")["Sample Name"]]

    logger.info(f"Found {len(samples)} samples. Outputting to {output}.")
    with open(output, "w") as f:
        for s in samples:
            f.write(s+'\n')


@click.command()
@click.argument("osd-id")
@click.option("--includes-assay-type", default=None)
@click.option("--includes-assay-type-on-platform", default=None)
@click.option("--includes-file-pattern", default=None)
@click.option("--excludes-file-pattern", default=None)
def check_if(osd_id, includes_assay_type, includes_assay_type_on_platform, includes_file_pattern, excludes_file_pattern):
    failed_criteria = None
    logger.info(f"Fetching information for {osd_id}")


    files = commons.get_table_of_files(osd_id)
    logger.info(f"Found {len(files)} files associated to {osd_id}")
    isa_path = isa.download_isa(accession = osd_id)

    try: # Try here is ensure isa path is deleted regardless of success
        ### Describe/Filter on Assay Types
        if includes_assay_type and failed_criteria == None:
            query_measurement, query_technology = includes_assay_type.split(",")
            assays_subtable = isa_archive.isa_investigation_subtables(Path(isa_path))["STUDY ASSAYS"]

            logger.trace(assays_subtable)
            # Add OSD ID
            assays_subtable["OSD_ID"] = osd_id

            # Add GLDS-ID
            # Assumes all files have GLDS number as prefix
            assays_subtable["GLDS_ID"] = files.iloc[0]['file_name'].split("_")[0]

            assays_subtable.set_index(["OSD_ID","GLDS_ID"], inplace = True)

            logger.trace(f"Found assay types: {assays_subtable}")

            # Check if desired assay type exists
            if assays_subtable.loc[
                            (assays_subtable["Study Assay Measurement Type"] == query_measurement) & 
                            (assays_subtable["Study Assay Technology Type"] == query_technology)
                            ].empty:
                failed_criteria = f"Could not find {query_measurement},{query_technology} in:\n {assays_subtable.to_dict(orient='records')}"
        
        ### Describe/Filter on Assay Types AND platform
        if includes_assay_type_on_platform and failed_criteria == None:
            query_measurement, query_technology, query_platform = includes_assay_type_on_platform.split(",")
            assays_subtable = isa_archive.isa_investigation_subtables(Path(isa_path))["STUDY ASSAYS"]

            logger.trace(assays_subtable)
            # Add OSD ID
            assays_subtable["OSD_ID"] = osd_id

            # Add GLDS-ID
            # Assumes all files have GLDS number as prefix
            assays_subtable["GLDS_ID"] = files.iloc[0]['file_name'].split("_")[0]

            assays_subtable.set_index(["OSD_ID","GLDS_ID"], inplace = True)

            logger.trace(f"Found assay types: {assays_subtable}")

            # Check if desired assay type exists
            if assays_subtable.loc[
                            (assays_subtable["Study Assay Measurement Type"].str.match(query_measurement)) & 
                            (assays_subtable["Study Assay Technology Type"].str.match(query_technology)) &
                            (assays_subtable["Study Assay Technology Platform"].str.match(query_platform))
                            ].empty:
                failed_criteria = f"Could not find {includes_assay_type_on_platform} in:\n {assays_subtable.to_dict(orient='records')}"


        ### Describe/Filter on files
        if includes_file_pattern and failed_criteria == None:
            logger.info(f"Searching for {len(files)} files")

            logger.trace(files)

            if len(commons.find_matching_filenames(
                            accession = osd_id, 
                            filename_pattern = includes_file_pattern
                            )
                    ) == 0:
                failed_criteria = (f"No files matching pattern '{includes_file_pattern}' could be located.")

        ### Describe/Filter on files
        if excludes_file_pattern and failed_criteria == None:
            logger.info(f"Searching for {len(files)} files")

            logger.trace(files)

            if len(commons.find_matching_filenames(
                            accession = osd_id, 
                            filename_pattern = excludes_file_pattern
                            )
                    ) != 0:
                failed_criteria = (f"Found files matching pattern '{excludes_file_pattern}' could be located.")

        # Check if any failed critera
        if failed_criteria is None:
            logger.success(f"{osd_id} matches supplied criteria!")
        else:
            logger.error(failed_criteria)
            sys.exit(-1)
    finally:
        # Teardown isa path regardless
        logger.info(f"Clean up: Removing {isa_path}")
        Path(isa_path).unlink()


osd.add_command(download_files)
osd.add_command(check_if)
osd.add_command(get_samples)

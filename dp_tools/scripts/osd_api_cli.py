from pathlib import Path
import sys

import click
from loguru import logger

from dp_tools.glds_api import commons, isa
from dp_tools.core.files import isa_archive

@click.group()
def osd():
    pass


@click.command()
@click.argument("osd-id")
@click.option("--includes-assay-type", default=None)
@click.option("--includes-file-pattern", default=None)
@click.option("--excludes-file-pattern", default=None)
def check_if(osd_id, includes_assay_type, includes_file_pattern, excludes_file_pattern):
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


osd.add_command(check_if)


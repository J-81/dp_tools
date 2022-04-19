import argparse
from pathlib import Path
import re
from typing import List
from dp_tools.bulkRNASeq.loaders import load_BulkRNASeq_STAGE_00
from dp_tools.core.configuration import load_full_config
from dp_tools.components.components import BulkRNASeqMetadataComponent
from dp_tools.glds_api.files import get_urls
import pandas as pd

import logging

log = logging.getLogger(__name__)

# TODO: refactor this with the analogous metadata component method
def isa_investigation_subtables(ISAarchive: Path) -> dict[str, pd.DataFrame]:
    tables: dict[str, pd.DataFrame] = dict()

    # track sub table lines
    table_lines: List[list] = list()
    key: str = None  # type: ignore

    try:
        [i_file] = (
            f
            for f in BulkRNASeqMetadataComponent.fetch_isa_files_external(ISAarchive)
            if f.name.startswith("i_")
        )
    except ValueError:
        raise FileNotFoundError(
            f"Could not find an i_* file inside: {ISAarchive.name}, is this an ISA archive?"
        )
    with open(i_file, "r") as f:
        for line in [l.rstrip() for l in f.readlines()]:
            # search for header
            if line in BulkRNASeqMetadataComponent._ISA_INVESTIGATION_HEADERS:
                if key != None:
                    tables[key] = pd.DataFrame(
                        table_lines
                    ).T  # each subtable is transposed in the i_file
                    table_lines = list()
                key = line  # set next table key
            else:
                tokens = line.split("\t")  # tab separated
                table_lines.append(tokens)
    tables[key] = pd.DataFrame(
        table_lines
    ).T  # each subtable is transposed in the i_file

    # reformat each table
    def clean_quotes(string: str) -> str:
        SINGLE_OR_DOUBLE_QUOTES = "\"'"
        # don't perform on non-string elements
        if not isinstance(string, str):
            return string
        else:
            return string.lstrip(SINGLE_OR_DOUBLE_QUOTES).rstrip(
                SINGLE_OR_DOUBLE_QUOTES
            )

    df: pd.DataFrame
    for key, df in tables.items():

        # note: as a ref, no reassign needed
        tables[key] = (
            df.rename(columns=df.iloc[0]).drop(df.index[0]).applymap(clean_quotes)
        )

    # ensure all expected subtables present
    assert set(tables.keys()) == BulkRNASeqMetadataComponent._ISA_INVESTIGATION_HEADERS

    return tables


def get_assay_table_path(
    ISAarchive: Path, configuration: dict, return_index: bool = False
) -> Path:
    """Retrieve the assay table file name that determined as a valid assay based on configuration.
    Specifically, defined in subsection 'ISA meta'

    :param study_assay_table: From the investigation file
    :type study_assay_table: pd.DataFrame
    :param configuration: Standard assay parsed config
    :type configuration: dict
    :return: Path to the found assay table
    :rtype: Path
    """
    config = configuration["ISA Meta"]
    # retrieve study assay subtable from I_file
    df = isa_investigation_subtables(ISAarchive)["STUDY ASSAYS"]

    # get valid tuples of measurement and technology types from configuration
    valid_measurements_and_technology_types: list[tuple[str, str]] = [
        (entry["measurement"], entry["technology"])
        for entry in config["Valid Study Assay Technology And Measurement Types"]
    ]

    # check for matching rows based on configuration tuple
    # one and only one row should match
    # not very efficient, but table should never be too large for this to be of concern
    matches: list[Path] = list()
    for valid_combination in valid_measurements_and_technology_types:
        log.debug(f"Searching subtable for {valid_combination}")
        match_row = df.loc[
            (
                df[["Study Assay Measurement Type", "Study Assay Technology Type"]]
                == valid_combination
            ).all(axis="columns")
        ]
        match_file = [Path(val) for val in match_row["Study Assay File Name"].values]
        matches.extend(match_file)

    # guard, one and only one should match
    assert (
        len(matches) == 1
    ), f"One and only one should match, instead got these matches: {matches}"

    # load assay table
    assay_file_path = matches[0]
    [assay_path] = [
        f
        for f in BulkRNASeqMetadataComponent.fetch_isa_files_external(ISAarchive)
        if f.name == assay_file_path.name
    ]

    if return_index:
        for valid_combination in valid_measurements_and_technology_types:
            [match_index] = df.loc[
                (
                    df[["Study Assay Measurement Type", "Study Assay Technology Type"]]
                    == valid_combination
                ).all(axis="columns")
            ].index.values
            return match_index

    return assay_path


def _parse_args():
    """Parse command line args."""
    parser = argparse.ArgumentParser(
        description=f"Script for downloading latest ISA from GLDS repository"
    )
    parser.add_argument(
        "--accession", metavar="GLDS-001", required=True, help="GLDS accession number"
    )
    parser.add_argument(
        "--config", metavar="0", default="Latest", help="Packaged config to use"
    )
    parser.add_argument(
        "--isa-archive",
        required=True,
        help="Local location of ISA archive file. Can be downloaded from the GLDS repository with 'dpt-get-isa-archive'",
    )

    args = parser.parse_args()
    return args


def main():
    args = _parse_args()
    isa_to_runsheet(args.accession, Path(args.isa_archive), str(args.config))


def isa_to_runsheet(accession: str, isa_archive: Path, config: str):
    configuration = load_full_config(config=config)
    i_tables = isa_investigation_subtables(isa_archive)
    a_table = pd.read_csv(
        get_assay_table_path(ISAarchive=isa_archive, configuration=configuration),
        sep="\t",
    )
    a_study_assays_index = get_assay_table_path(
        ISAarchive=isa_archive, configuration=configuration, return_index=True
    )
    [s_file] = (
        f
        for f in BulkRNASeqMetadataComponent.fetch_isa_files_external(isa_archive)
        if f.name.startswith("s_")
    )
    s_table = pd.read_csv(s_file, sep="\t")
    df_merged = s_table.merge(a_table, on="Sample Name").set_index(
        "Sample Name", drop=True
    )
    df_final = pd.DataFrame(index=df_merged.index)
    # extract from Investigation table first
    investigation_source_entries = [
        entry
        for entry in configuration["Staging"]["General"]["Required Metadata"][
            "From ISA"
        ]
        if entry["ISA Table Source"] == "Investigation"
    ]
    for entry in investigation_source_entries:
        # handle special cases
        if entry.get("True If Includes At Least One"):
            overlap = set(entry["True If Includes At Least One"]).intersection(
                set(i_tables[entry["Investigation Subtable"]][entry["ISA Field Name"]])
            )
            df_final[entry["Runsheet Column Name"]] = bool(overlap)
            continue

        target_investigation_column = i_tables[entry["Investigation Subtable"]].loc[
            a_study_assays_index
        ]
        df_final[entry["Runsheet Column Name"]] = target_investigation_column[
            entry["ISA Field Name"]
        ]

    # extract from assay table first
    assay_source_entries = [
        entry
        for entry in configuration["Staging"]["General"]["Required Metadata"][
            "From ISA"
        ]
        if entry["ISA Table Source"] in ["Assay", "Sample", ["Assay", "Sample"]]
        and entry.get("Autoload", True) != False
    ]
    for entry in assay_source_entries:
        assert list(df_final.index) == list(df_merged.index)
        if entry.get("Runsheet Index"):
            # already set and checked above
            continue
        else:
            # factor value style extraction
            if entry.get("Matches Multiple Columns") and entry.get("Match Regex"):
                # find matching columns
                match_cols = [
                    (i, col, df_merged[col])
                    for i, col in enumerate(df_merged.columns)
                    if re.match(pattern=entry.get("Match Regex"), string=col)
                ]

                # check if columns require appending unit
                if entry.get("Append Column Following"):
                    match_i: int  # index in matching column list
                    df_i: int  # index in merged dataframe
                    col: str
                    original_series: pd.Series
                    for match_i, (df_i, col, original_series) in enumerate(match_cols):
                        # scan through following columns
                        for scan_col in df_merged.iloc[:, df_i:].columns:
                            # check if another 'owner' column is scanned, this means Unit was not found
                            if any(
                                [
                                    scan_col.startswith("Parameter Value["),
                                    scan_col.startswith("Factor Value["),
                                    scan_col.startswith("Characteristics ["),
                                ]
                            ):
                                break
                            if scan_col == entry.get("Append Column Following"):
                                resolved_series = original_series + df_merged[scan_col]
                                match_cols[match_i] = df_i, col, resolved_series
                                break

                # finally add this information into dataframe
                for _, col_name, series in match_cols:
                    df_final[col_name] = series
            else:
                # CAUTION: normally this wouldn't be safe as the order of rows isn't considered.
                # In this block, the indices are checked for parity already making this okay
                match entry["ISA Field Name"]:
                    case str():
                        series_to_add = df_merged[entry["ISA Field Name"]]
                    # handles cases where the field name varies
                    case list():
                        [found_col] = (
                            col
                            for col in df_merged.columns
                            if col in entry["ISA Field Name"]
                        )
                        series_to_add = df_merged[found_col]

                if entry.get("Remapping"):
                    df_final[entry["Runsheet Column Name"]] = series_to_add.map(
                        lambda val: entry.get("Remapping")[val]
                    )
                else:
                    df_final[entry["Runsheet Column Name"]] = series_to_add

    # now find glds urls and append
    user_source_entries = [
        entry
        for entry in configuration["Staging"]["General"]["Required Metadata"][
            "From User"
        ]
        if entry.get("Data Asset Keys")  # indicates an exact name can be determined
    ]

    urls = get_urls(accession=accession)
    for entry in user_source_entries:
        putative_filenames: set = set()
        # iterate through keys
        for key in entry["Data Asset Keys"]:
            data_asset_conf = configuration["data assets"][key]
            template = (
                configuration["ISA Meta"]["Global file prefix"].format(
                    datasystem=accession
                )
                + data_asset_conf["processed location"][-1]
            )
            putative_filenames = putative_filenames.union(
                set(
                    df_final.index.map(
                        lambda sample_name: template.format(sample=sample_name)
                    )
                )
            )

            # check if these are listed
            mapped_to_url = list(
                map(
                    lambda fn: urls.get(fn, dict()).get("url", False),
                    putative_filenames,
                )
            )
            if all(mapped_to_url):
                df_final[entry["Runsheet Column Name"]] = mapped_to_url

    # output file
    output_fn = (
        f"{accession}_{configuration['NAME']}_v{configuration['VERSION']}_runsheet.csv"
    )
    log.info(f"Writing runsheet to: {output_fn} with {len(df_final)} columns")
    df_final.to_csv(output_fn)

    return df_final


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()

""" Functions that take a data model and produce reporter files 
Depends on standard data asset metadata as loaded from packaged config files.
"""
from collections import defaultdict
import os
from pathlib import Path
import re
from typing import Union
from dp_tools.core.entity_model import TemplateDataset
import pkg_resources


import logging

log = logging.getLogger(__name__)

from schema import Schema, Or, And
import yaml
import pandas as pd


def load_config(
    config: Union[tuple[str, str], Path], subsection: str = "ISA Meta"
) -> dict:
    if isinstance(config, tuple):
        configuration = yaml.safe_load(
            pkg_resources.resource_string(
                __name__,
                os.path.join("..", "config", f"{config[0]}_v{config[1]}.yaml"),
            )
        )
    elif isinstance(config, Path):
        configuration = yaml.safe_load(config.open())

    # filter to relevant subsection
    sub_configuration = configuration[subsection]

    log.debug("Loaded the following validation config: {conf_validation}")

    # validate with schema
    config_schema = Schema(
        {
            "Valid Study Assay Technology And Measurement Types": [
                {"measurement": str, "technology": str}
            ],
            "Global file prefix": str,
        }
    )

    config_schema.validate(sub_configuration)

    return sub_configuration


def unmangle_columns(columns: list[str]) -> list[str]:
    """Utility function to convert "X.1...X.N" into "X...X", reversing the normal column name mangle
    ref: https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html

    :param columns: Original mangled column names
    :type columns: list[str]
    :return: Unmangled columns names in same order
    :rtype: list[str]
    """

    def remove_mangle(s: str) -> str:
        _RE_MANGLE = r"\.\d+$"  # e.g. {".1",".18"}
        return re.sub(pattern=_RE_MANGLE, repl="", string=s)

    new_cols = [remove_mangle(col) for col in columns]
    return new_cols


def update_curation_tables(
    dataset: TemplateDataset,
    config: Union[tuple[str, str], Path],
    output_file: str = None,
):
    """Updates existing curation investigation and assay tables with processed data

    :param dataset: A loaded dataset object
    :type dataset: TemplateDataset
    :param config: A assay type configuration file specifier formatted as follows: ('assay','version') or local path to configuration file
    :type config: Union[tuple[str, str], Path]
    :param output_file: The name of the updated output tables, defaults to using the same name as the original tables
    :type output_file: str, optional
    """
    _DATA_ASSET_COL_PREFIX = "Parameter Value["
    _DATA_ASSET_COL_SUFFIX = "]"

    configuration = load_config(config)

    # retrieve study assay subtable from I_file
    df = dataset.metadata.isa_investigation_subtables["STUDY ASSAYS"]

    # get valid tuples of measurement and technology types from configuration
    valid_measurements_and_technology_types: list[tuple[str, str]] = [
        (entry["measurement"], entry["technology"])
        for entry in configuration["Valid Study Assay Technology And Measurement Types"]
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
    assay_file_name = matches[0]
    [assay_path] = [
        f for f in dataset.metadata.fetch_isa_files() if f.name == assay_file_name.name
    ]

    df_assay = pd.read_csv(assay_path, sep="\t").set_index(keys="Sample Name")

    # e.g.
    # {"Parameter Value[New1/subnew1]" : {"sample1":"asset",...}}

    # generate table of sample wise files
    # associating all non-samplewise files to ALL samples
    all_samples = list(dataset.samples.keys())

    new_cols: dict[str, dict[str, str],] = defaultdict(
        lambda: defaultdict(set)  # type: ignore
    )

    column_order: dict = dict()

    for asset in dataset.all_data_assets:
        if not asset.metadata["resource categories"]["publish to repo"]:
            continue

        # alias
        resource_config = asset.metadata["resource categories"]

        # format header
        category_string = (
            f"{resource_config['subcategory']}/{resource_config['subdirectory']}"
            if all(
                [
                    resource_config["subdirectory"] != "",
                    resource_config["include subdirectory in table"],
                ]
            )
            else resource_config["subcategory"]
        )
        header = _DATA_ASSET_COL_PREFIX + category_string + _DATA_ASSET_COL_SUFFIX
        associated_samples = (
            [asset.metadata["template_kwargs"].get("sample", None)]
            if asset.metadata["template_kwargs"].get("sample", None)
            else all_samples
        )

        # TODO: this likely will be better placed elsewhere but for now this is fine
        # Track column order here as a dict
        column_order[header] = asset.metadata["resource categories"]["table order"]

        # associate file names to each sample
        for sample in associated_samples:
            # format file names
            table_filename = (
                configuration["Global file prefix"].format(
                    datasystem=dataset.dataSystem.name
                )
                + asset.path.name
            )

            # let the entries start as lists
            # these will be converted to sorted and joined strings at the end
            new_cols[header][sample].add(table_filename)  # type: ignore

    # convert all sets by
    #  conversion to list
    #  sorting alphabetically
    #  joining by comma
    for header, header_wise in new_cols.items():
        for sample, sample_wise in header_wise.items():
            new_value = ", ".join(sorted(list(new_cols[header][sample])))

            new_cols[header][sample] = new_value

    # TODO:
    # sort columns based on configuration
    # check that all column order values are unique and non-negative (the value for unpublished items)
    order_values: set[int] = set()
    for col, order in column_order.items():
        assert order >= 0, f"Order value must be non-negative: ({order})"
        assert (
            order not in order_values
        ), f"Duplicate column order value ({order}) found for header: {col}"
        order_values.add(order)

    # original columns in order
    orig_columns = list(df_assay.columns)
    sorted_new_columns = [
        col for col in sorted(column_order, key=lambda k: column_order[k])
    ]

    df_assay_extended = df_assay.join(pd.DataFrame(new_cols))

    # now reorder with both new and original columns
    df_assay_extended = df_assay_extended[orig_columns + sorted_new_columns]

    # guards
    assert len(df_assay_extended.index) == len(
        df_assay.index
    ), f"After join, index length did not stay the same: old_length-{len(df_assay.index)} new_length-{len(df_assay_extended.index)}"
    assert len(df_assay_extended.columns) > len(
        df_assay.columns
    ), f"After join, no new columns were added"

    # create default output file name if not provided
    output_file = (
        output_file
        if output_file is not None
        else f"{dataset.name}_assay_table_after_processing.tsv"
    )

    # dropped NA columns
    # these often exist due to an old requirement to include ontology columns
    # even when no such ontology values were present
    df_assay_extended = df_assay_extended.dropna(axis="columns", how="all")

    # used unmangled columns aliases when writing to file
    df_assay_extended.to_csv(
        output_file, header=unmangle_columns(df_assay_extended.columns), sep="\t"
    )

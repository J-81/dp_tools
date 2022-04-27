""" Functions that take a data model and produce reporter files 
Depends on standard data asset metadata as loaded from packaged config files.
"""
from collections import defaultdict
import os
from pathlib import Path
import re
from typing import Optional, Union
from dp_tools.core.entity_model import TemplateDataset
import pkg_resources


import logging

log = logging.getLogger(__name__)

from schema import Schema, Or, And
import yaml
import pandas as pd

# constants
_PARAMETER_VALUE_COL_PREFIX = "Parameter Value["
_PARAMETER_VALUE_COL_SUFFIX = "]"
_PROTOCOL_REP_COL_PREFIX = "Protocol REF"


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
            "Post Processing Add Study Protocol": str,
        }
    )

    config_schema.validate(sub_configuration)

    return sub_configuration


def load_ISA_investigation_config() -> dict:
    configuration = yaml.safe_load(
        pkg_resources.resource_string(
            __name__,
            os.path.join("..", "config", f"ISA_investigation.yaml"),
        )
    )

    log.debug("Loaded the ISA investigation config")

    # TODO: validate with schema
    # config_schema.validate(sub_configuration)

    return configuration


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


def get_assay_table_path(dataset: TemplateDataset, configuration: dict) -> Path:
    """Retrieve the assay table file name that determined as a valid assay based on configuration.
    Specifically, defined in subsection 'ISA meta'

    :param dataset: A dataset object including a metadata component with an attached ISA archive data asset
    :type dataset: TemplateDataset
    :param configuration: Standard assay parsed config
    :type configuration: dict
    :return: Path to the found assay table
    :rtype: Path
    """
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
    assay_file_path = matches[0]
    [assay_path] = [
        f for f in dataset.metadata.fetch_isa_files() if f.name == assay_file_path.name
    ]

    return assay_path


def generate_new_column_dicts(
    dataset: TemplateDataset, configuration: dict
) -> tuple[dict, dict]:
    """Based on data assets in the dataset and configuration, generates dictionaries for dataframe extension.
    Specifically, the "data assets" configuration subsection is used.

    :param dataset: A dataset object
    :type dataset: TemplateDataset
    :param configuration: Standard assay parsed config
    :type configuration: dict
    :return: Two dictionaries, the first with dataframe ready new data assets and the second denoting the intended order of those new columns
    :rtype: tuple[dict, dict]
    """
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
        header = (
            _PARAMETER_VALUE_COL_PREFIX + category_string + _PARAMETER_VALUE_COL_SUFFIX
        )
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

    order_values: set[int] = set()
    for col, order in column_order.items():
        assert order >= 0, f"Order value must be non-negative: ({order})"
        assert (
            order not in order_values
        ), f"Duplicate column order value ({order}) found for header: {col}"
        order_values.add(order)

    return new_cols, column_order


def extend_assay_dataframe(
    df_orignal: pd.DataFrame, new_column_data: dict, new_column_order: dict
):
    # original columns in order
    orig_columns = list(df_orignal.columns)
    sorted_new_columns = [
        col for col in sorted(new_column_order, key=lambda k: new_column_order[k])
    ]
    df_extended = df_orignal.join(pd.DataFrame(new_column_data))

    # now reorder with both new and original columns
    df_extended = df_extended[orig_columns + sorted_new_columns]

    # guards
    assert len(df_extended.index) == len(
        df_orignal.index
    ), f"After join, index length did not stay the same: old_length-{len(df_orignal.index)} new_length-{len(df_extended.index)}"
    assert len(df_extended.columns) > len(
        df_orignal.columns
    ), f"After join, no new columns were added"

    # dropped NA columns
    # these often exist due to an old requirement to include ontology columns
    # even when no such ontology values were present
    df_extended = df_extended.dropna(axis="columns", how="all")

    return df_extended


def sync_investigation_and_assay_dataframes(
    df_investigation, column_contents, column_order, configuration
):
    configuration_ISA_investigation = load_ISA_investigation_config()

    print(1)


def get_parameter_values(
    df: pd.DataFrame, drop_ontology: bool = True
) -> dict[str, list]:
    """Extract the value from parameter value columns.
    Such that 'Parameter Value[someVal]' gives 'someVal'
    These are associated with the first preceeding
    "Protocol REF" column (which is likely mangled by pandas)

    :param df: A dataframe representing an assay or sample table from the ISA spec
    :type df: pd.DataFrame
    :param drop_ontology: In legacy data,
    :type drop_ontology: pd.DataFrame
    :return: A mapping of parameter values to their protocol
    :rtype: dict[str, list]
    """

    def _extract_value(s: str) -> str:
        return s.removeprefix(_PARAMETER_VALUE_COL_PREFIX).removesuffix(
            _PARAMETER_VALUE_COL_SUFFIX
        )

    result = defaultdict(list)

    for col in df.columns:
        if col.startswith(_PROTOCOL_REP_COL_PREFIX):
            [cur_protocol] = df[col].unique()
        if col.startswith(_PARAMETER_VALUE_COL_PREFIX):
            result[cur_protocol].append(_extract_value(col))

    return result


def add_protocol(
    df_investigation: pd.DataFrame,
    protocol_key: str,
    df_assay: Optional[pd.DataFrame] = None,
):
    # retrieve study protocol
    conf_isa = load_ISA_investigation_config()
    target_protocol = conf_isa["STUDY PROTOCOLS"][protocol_key]

    # get parameter values
    if df_assay is not None:
        params_assay = get_parameter_values(df_assay, drop_ontology=True)
        # compare against target protocol expected 'Study Protocol Parameters Name' list
        print(1)
    print(1)


def setup_output_target(
    output_file: Optional[str], original_path: Path, output_dir: str = "updated_curation_tables"
) -> Path:
    """Set ups target output file location.  Uses specified output_file name if provided.
    Defaults to the original table filename but saves to a directory 'output_dir' in either case.

    :param output_file: Specifies an alternative output filename for the extended table
    :type output_file: Optional[str]
    :param original_path: Specifies the orignal path, used to determine a default filename.
    :type original_path: Path
    :param output_dir: Specifies the name of the directory for the output file, defaults to "updated_tables"
    :type output_dir: str, optional
    :return: The output target location.
    :rtype: Path
    """
    # create default output file name if not provided

    final_output_target: Path = (
        Path(output_file)
        if output_file is not None
        else Path(output_dir) / original_path.name
    )
    Path(output_dir).mkdir(exist_ok=True)

    return final_output_target


def update_curation_tables(
    dataset: TemplateDataset,
    config: Union[tuple[str, str], Path],
    output_file: str = None,
    investigation_table: bool = False,
):
    """Updates existing curation investigation and assay tables with processed data.
    This function extends the existing assay table with processing data assets.
    The naming and order of the new columns is configured in the data assets config section.

    :param dataset: A loaded dataset object
    :type dataset: TemplateDataset
    :param config: A assay type configuration file specifier formatted as follows: ('assay','version') or local path to configuration file
    :type config: Union[tuple[str, str], Path]
    :param output_file: The name of the updated output tables, defaults to using the same name as the original tables
    :type output_file: str, optional
    :param investigation_table: Controls the updating and syncing of the investigation table, defaults to False
    :type investigation_table: bool, optional
    """
    configuration = load_config(config)

    # first extend the assay table
    # load assay dataframe
    assay_path = get_assay_table_path(dataset, configuration)
    df_assay = pd.read_csv(assay_path, sep="\t").set_index(keys="Sample Name")

    column_contents, column_order = generate_new_column_dicts(dataset, configuration)

    df_assay_extended = extend_assay_dataframe(df_assay, column_contents, column_order)

    # now modify investigation table accordingly
    assay_path = get_assay_table_path(dataset, configuration)

    if investigation_table:
        raise NotImplementedError(
            "Updating Investigation table is not currently implemented."
        )
        df_investigation = dataset.metadata.isa_investigation_subtables[
            "STUDY PROTOCOLS"
        ]
        resolved_protocol_key = configuration[
            "Post Processing Add Study Protocol"
        ].format(LEADCAP_organism=dataset.metadata.organism.capitalize())
        df_investigation = add_protocol(
            df_investigation, protocol_key=resolved_protocol_key, df_assay=df_assay
        )
        df_investigation_extended = sync_investigation_and_assay_dataframes(
            df_investigation, column_contents, column_order, configuration
        )
        # NOT Implemented: write_investigation_table()

    final_output_target = setup_output_target(output_file, original_path=assay_path)

    # used unmangled columns aliases when writing to file
    df_assay_extended.to_csv(
        final_output_target,
        header=unmangle_columns(df_assay_extended.columns),
        sep="\t",
    )

""" Functions that take a data model and produce reporter files 
Depends on standard data asset metadata as loaded from packaged config files.
"""
from collections import defaultdict
import hashlib
import os
from pathlib import Path
import re
from typing import Optional, TypedDict, Union
from dp_tools.core.configuration import load_config
import pkg_resources


import logging

from dp_tools.core.entity_model import Dataset, Group, Sample
from dp_tools.core.files import isa_archive

log = logging.getLogger(__name__)

from schema import Schema
import yaml
import pandas as pd

# constants
_PARAMETER_VALUE_COL_PREFIX = "Parameter Value["
_PARAMETER_VALUE_COL_SUFFIX = "]"
_PROTOCOL_REP_COL_PREFIX = "Protocol REF"


def _load_config(
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
    if subsection == "ISA Meta":
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


def get_assay_table_path(isaArchive: Path, configuration: dict) -> Path:
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
    df = isa_archive.isa_investigation_subtables(isaArchive)["STUDY ASSAYS"]

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
        f
        for f in isa_archive.fetch_isa_files(isaArchive)
        if f.name == assay_file_path.name
    ]

    return assay_path


# typehints mapping to config
ResourceCategory = TypedDict(
    "ResourceCategory",
    {
        "subcategory": "str",
        "subdirectory": "str",
        "publish to repo": "bool",
        "include subdirectory in table": "bool",
        "table order": "int",
    },
)

DataAssetConfig = TypedDict(
    "DataAssetConfig",
    {"processed location": "list", "resource categories": "ResourceCategory"},
)


def get_repolike_category_string(data_asset_metadata: ResourceCategory) -> str:
    if all(
        [
            data_asset_metadata["subdirectory"] != "",
            data_asset_metadata["include subdirectory in table"],
        ]
    ):
        return f"{data_asset_metadata['subcategory']}/{data_asset_metadata['subdirectory']}"
    else:
        return data_asset_metadata["subcategory"]


def generate_new_column_dicts(
    dataset: Dataset, configuration: dict
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

    new_cols: dict[str, dict[str, str],] = defaultdict(
        lambda: defaultdict(set)  # type: ignore
    )

    column_order: dict = dict()

    for asset in dataset.get_assets():
        if not asset.config["resource categories"]["publish to repo"]:
            continue

        # format header
        category_string = get_repolike_category_string(
            asset.config["resource categories"]
        )
        header = (
            _PARAMETER_VALUE_COL_PREFIX + category_string + _PARAMETER_VALUE_COL_SUFFIX
        )

        # sample names here are based on processing sample names!
        match asset.owner:
            case Sample():
                associated_samples = [asset.owner.name]  # Must be a list
            case Dataset():
                associated_samples = list(dataset.samples.keys())
            case Group():
                raise NotImplementedError(
                    "Group ownership is not implemented in current data model"
                )

        # now remap those processing sample names to their orignal names,
        # required for fing to orignal assay table
        processing_to_orignal_mapping = pd.read_csv(
            dataset.data_assets["runsheet"].path, index_col="Sample Name"
        )["Original Sample Name"].to_dict()
        # TODO: ineffecient, should only iterate once
        remapped_samples = [
            f"{sample} remapped to {processing_to_orignal_mapping[sample]}"
            for sample in associated_samples
            if sample != processing_to_orignal_mapping[sample]
        ]
        associated_samples = [
            processing_to_orignal_mapping[sample] for sample in associated_samples
        ]
        if remapped_samples:
            log.info(
                f"Post processing using remapped samples for the following: {remapped_samples} for header: '{header}'"
            )

        # TODO: this likely will be better placed elsewhere but for now this is fine
        # Track column order here as a dict
        column_order[header] = asset.config["resource categories"]["table order"]

        # associate file names to each sample
        for sample in associated_samples:
            # format file names
            table_filename = (
                configuration["Global file prefix"].format(datasystem=dataset.name)
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
            new_value = ",".join(sorted(list(new_cols[header][sample])))

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
    output_file: Optional[str],
    original_path: Path,
    output_dir: str = "updated_curation_tables",
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
    dataset: Dataset,
    config: Union[tuple[str, str], Path],
    output_file: str = None,
    investigation_table: bool = False,
) -> pd.DataFrame:
    """Updates existing curation investigation and assay tables with processed data.
    This function extends the existing assay table with processing data assets.
    The naming and order of the new columns is configured in the data assets config section.

    :param dataset: A loaded dataset object
    :type dataset: Dataset
    :param config: A assay type configuration file specifier formatted as follows: ('assay','version') or local path to configuration file
    :type config: Union[tuple[str, str], Path]
    :param output_file: The name of the updated output tables, defaults to using the same name as the original tables
    :type output_file: str, optional
    :param investigation_table: Controls the updating and syncing of the investigation table, defaults to False
    :type investigation_table: bool, optional
    """
    configuration = _load_config(config)

    # first extend the assay table
    # load assay dataframe
    isaArchive = dataset.data_assets["ISA Archive"].path
    assay_path = get_assay_table_path(isaArchive, configuration)
    df_assay = pd.read_csv(assay_path, sep="\t").set_index(keys="Sample Name")

    column_contents, column_order = generate_new_column_dicts(dataset, configuration)

    df_assay_extended = extend_assay_dataframe(df_assay, column_contents, column_order)

    # now modify investigation table accordingly
    assay_path = get_assay_table_path(isaArchive, configuration)

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

    return df_assay_extended


class Md5sum_row(TypedDict):
    resource_category: str
    filename: str
    md5sum: str


class Md5sum_row_with_tags(Md5sum_row):
    tags: str


def compute_md5sum(file_path: Path) -> str:
    return hashlib.md5(file_path.open("rb").read()).hexdigest()


# ALLOWED MISSING KEYS CONSTANTS
ALLOWED_MISSING_KEYS_FOR_SINGLE_END = {
    "raw forward reads fastq GZ",
    "raw reverse reads fastq GZ",
    "forward reads trimming report",
    "reverse reads trimming report",
    "trimmed reverse reads fastq GZ",
    "trimmed forward reads fastq GZ",
    "inner distance MultiQC directory ZIP",
}
""" Data assets unique to paired end layouts """

ALLOWED_MISSING_KEYS_FOR_PAIRED_END = {
    "raw reads fastq GZ",
    "trimmed reads fastq GZ",
    "reads trimming report",
}
""" Data assets unique to single end layouts """


ALLOWED_MISSING_KEYS_FOR_NON_ERCC = {
    "ERCC normalized DESeq2 contrasts table",
    "ERCC normalized DESeq2 normalized counts table",
    "ERCC normalized DESeq2 annotated DGE table",
    "ERCC sample table",
    "ERCC analysis HTML",
}
""" Data assets unique to ERCC spiked-in datasets """


def generate_md5sum_table(
    dataset: Dataset,
    config: Union[tuple[str, str], Path],
    allowed_unused_keys: set[str] = None,
    include_tags: bool = False,
) -> pd.DataFrame:
    # set default empty set
    if allowed_unused_keys is None:
        allowed_unused_keys = set()

    loaded_config = load_config(config)

    # generate all md5sums for all data assets that will be published
    # table columns: resource_category, filename, md5sum
    data: Union[list[Md5sum_row_with_tags], list[Md5sum_row]] = list()
    for asset in dataset.get_assets():
        if not asset.config["resource categories"]["publish to repo"]:
            continue

        # catch rare cases where a data asset is 'psuedo loaded'
        # this occurs when a data asset is loaded without verifying
        # the asset exists. Used when a data asset is part of the repo publish assets
        # but not generated during automated processing and rather added after the fact
        if asset.putative:
            data.append(
                {
                    "resource_category": get_repolike_category_string(
                        asset.config["resource categories"]
                    ),
                    "filename": asset.path.name,
                    "md5sum": "USER MUST ADD MANUALLY!",
                }
                | ({"tags": asset.config["tags"]} if include_tags else {})
            )
            continue  # to next data asset

        # branch for data asset files
        if asset.path.is_file():
            data.append(
                {
                    "resource_category": get_repolike_category_string(
                        asset.config["resource categories"]
                    ),
                    "filename": asset.path.name,
                    "md5sum": compute_md5sum(asset.path),
                }
                | ({"tags": asset.config["tags"]} if include_tags else {})
            )
        # branch for data asset dirs
        elif asset.path.is_dir():
            for sub_asset in asset.path.iterdir():
                data.append(
                    {
                        "resource_category": get_repolike_category_string(
                            asset.config["resource categories"]
                        ),
                        "filename": sub_asset.name,
                        "md5sum": compute_md5sum(sub_asset),
                    }
                    | ({"tags": asset.config["tags"]} if include_tags else {})
                )

    # compare all data assets against configuration
    asset_keys_in_dataset: set[str] = {asset.key for asset in dataset.get_assets()}
    publishable_asset_keys_in_config: set[str] = {
        key
        for key, value in loaded_config["data assets"].items()
        if value["resource categories"]["publish to repo"]
    }

    missing_publishables_by_key = (
        publishable_asset_keys_in_config - asset_keys_in_dataset
    )

    # report any missing data assets that are supposed to be published
    if not_allowed_missing := missing_publishables_by_key - allowed_unused_keys:
        raise ValueError(
            f"The following data assets keys are not allowed to be missing: {not_allowed_missing}. Allowed missing keys: {allowed_unused_keys}."
        )

    # generate dataframe and return
    df = (
        pd.DataFrame(data)
        .sort_values(by=["resource_category", "filename"])
        .drop_duplicates(subset=["md5sum", "resource_category", "filename"])
    )
    return df

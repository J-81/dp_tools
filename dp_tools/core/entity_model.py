from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from string import Formatter
import uuid
from typing import TypedDict, Union
import logging

from dp_tools.core.check_model import FlagCode

log = logging.getLogger(__name__)

import pandas as pd
import multiqc

# MULTIQC MONKEY PATCH TO ADDRESS ISSUE: https://github.com/ewels/MultiQC/issues/1643
multiqc.config.logger.hasHandlers = (
    lambda: False
)  # this means the logger never gets purged, but more importantly prevents a log purge based exceptoin


from dp_tools.core.utilites.multiqc_tools import (
    format_plots_as_dataframe,
    get_general_stats,
)


def get_id():
    uuid4 = uuid.uuid4()
    return str(uuid4)


#########################################################################
# DATAFILE
#########################################################################
ExperimentalEntity = Union["Dataset", "Group", "Sample"]


@dataclass
class DataAsset:
    key: str
    """ Configuration key for this data asset"""
    path: Path
    """ Path object for the data asset """
    config: dict
    """ Configuration dict directly from yaml file """
    owner: ExperimentalEntity
    """ The owner of the data asset, an experimental entity """
    putative: bool = field(default=False)
    """ Indicates if the data asset is loaded putatively (i.e. expected to exist in the future"""


DataAssetDict = dict[str, DataAsset]


def dataSystem_from_runsheet(runsheet_path: Path) -> "DataSystem":
    # Parse runsheet name
    runsheet_meta = DataSystem.parse_runsheet_name(runsheet_path.name)

    return DataSystem(name=runsheet_meta["dataset_name"])


@dataclass
class DataSystem:
    name: str
    """Name of the DataSystem. e.g. 'GLDS-48'"""
    datasets: dict[str, "Dataset"] = field(default_factory=dict, repr=False)
    """Datasets objects indexed by their names"""

    class RunsheetMeta(TypedDict):
        dataset_name: str
        dataset_type: str
        version: str

    @property
    def dataset(self):
        assert (
            len(self.datasets) == 1
        ), "Can only call 'dataset' if only one dataset exists in the dataSystem"
        return list(self.datasets.values())[0]

    @staticmethod
    def parse_runsheet_name(runsheet_name: str) -> RunsheetMeta:
        """Parse runsheet name for runsheet metadata.

        :param runsheet_name: Runsheet file name
        :type runsheet_name: str
        :return: Dictionary containing data about the runsheet
        :rtype: RunsheetMeta
        """
        # remove '_runsheet.csv' suffix
        s = runsheet_name.replace("_runsheet.csv", "")

        # split by '_' and peel off last two tokens (this allows dataset names to include underscores)
        tokens = s.split("_")

        return {
            "dataset_name": "_".join(tokens[:-2]),
            "dataset_type": tokens[-2],
            "version": tokens[-1],
        }

    def dataset_from_runsheet(self, runsheet_path: Path) -> "Dataset":
        """Initiates a dataset  from runsheet and attaches this data system

        :param runsheet_path: Path to runsheet
        :type runsheet_path: Path
        """
        log.info(f"Loading dataset from {runsheet_path}")
        df = pd.read_csv(runsheet_path)

        # Parse runsheet name
        runsheet_meta = DataSystem.parse_runsheet_name(runsheet_path.name)

        # initate dataset
        dataset = Dataset(
            name=runsheet_meta["dataset_name"], type=runsheet_meta["dataset_type"]
        )

        # load samples in dataset
        for sample in df["Sample Name"]:
            dataset.samples[sample] = Sample(name=sample)

        # attach dataset to DataSystem
        self.datasets[dataset.name] = dataset

        # attach runsheet as dataset asset
        dataset.data_assets["runsheet"] = DataAsset(
            key="runsheet",
            path=runsheet_path,
            config={"processed location": runsheet_path},
            owner=dataset,
        )

        # Add metadata from runsheet_path
        dataset.metadata.update(df.loc[:, df.nunique() == 1].iloc[0].to_dict())

        # returns attach dataset for convienience
        return dataset


@dataclass
class Dataset:
    name: str
    type: str
    metadata: dict[str, str] = field(default_factory=dict)
    samples: dict[str, "Sample"] = field(default_factory=dict, repr=False)
    groups: dict[str, "Group"] = field(default_factory=dict, repr=False)
    data_assets: DataAssetDict = field(default_factory=dict, repr=False)
    ALLOWED_FORMAT_KEYS: tuple[str, str, str] = field(
        default=("dataset", "sample", "group"), repr=False
    )
    loaded_assets_dicts: list[dict] = field(default_factory=list, repr=False)

    def _load_asset(
        self,
        asset: Path,
        key: str,
        config: dict,
        owner: ExperimentalEntity,
        putative: bool,
    ) -> DataAsset:
        if "*" in asset.name:
            try:
                [asset] = asset.parent.glob(asset.name)
            except ValueError as exc:
                raise ValueError(
                    f"Failed to locate data asset using glob pattern: '{asset.name}'"
                ) from exc
        if not putative:
            assert asset.exists(), f"Failed to load asset at path '{asset}'"
            self.loaded_assets_dicts.append(
                {
                    "index": (owner.name, "Data Assets", key),
                    "description": f"Check data asset for key '{key}' exists",
                    "function": self._load_asset.__name__,
                    "code": FlagCode.GREEN,
                    "message": f"Data asset located: {asset.name}",
                    "code_level": FlagCode.GREEN.value,
                    "kwargs": {"asset": asset, "config": config},
                    "config": {},
                }
            )
        else:  # Confer as a RED flag for log purposes
            self.loaded_assets_dicts.append(
                {
                    "index": (owner.name, "Data Assets", key),
                    "description": f"Putative load for data asset for key '{key}' (This means the data asset is expected to exist in the future)",
                    "function": self._load_asset.__name__,
                    "code": FlagCode.RED,
                    "message": f"Future Data asset to be located: {asset.name}",
                    "code_level": FlagCode.RED.value,
                    "kwargs": {"asset": asset, "config": config},
                    "config": {},
                }
            )
        return DataAsset(
            key=key, path=asset, config=config, owner=owner, putative=putative
        )

    @property
    def loaded_assets_report(self) -> pd.DataFrame:
        return pd.DataFrame(self.loaded_assets_dicts).set_index(keys="index")

    # TODO: dict -> better typehint via typeddict
    def load_data_asset(
        self, data_asset_config: dict, root_dir: Path, name: str, putative: bool = False
    ):
        # Check if a dataset conditional preempts loading
        if conditions := data_asset_config.get("conditional on dataset", None):
            for condition in conditions:
                for metadata_key, allowed_values in condition.items():
                    if self.metadata[metadata_key] not in allowed_values:
                        log.warning(
                            f"Skipping loading of requested asset with key '{name}' due to failed dataset metadata condition: {conditions}"
                        )
                        return

        # Template Path
        location_template = Path(root_dir, *data_asset_config["processed location"])

        # Infer ownership level (e.g. dataset, group, sample)
        # Most fine grained scope is the owner, no template defaults to dataset
        format_keys = tuple(
            [
                tup[1]
                for tup in Formatter().parse(location_template.name)
                if tup[1] is not None
            ]
        )

        assert set(format_keys).issubset(
            self.ALLOWED_FORMAT_KEYS
        ), f"Found these template arguments: {format_keys} but all template arguments must be in: {self.ALLOWED_FORMAT_KEYS}"

        match format_keys:
            case ["sample"]:
                owner = "sample"
            case ["group"]:
                owner = "group"
            case ["dataset"]:
                owner = "dataset"
            case _:
                owner = "dataset"

        # Locate data asset
        match owner:
            case "dataset":
                unloaded_asset = Path(str(location_template).format(dataset=self.name))
                asset = self._load_asset(
                    unloaded_asset,
                    key=name,
                    config=data_asset_config,
                    owner=self,
                    putative=putative,
                )
                self.data_assets[name] = asset
            case "group":
                for group in self.groups.values():
                    unloaded_asset = Path(
                        str(location_template).format(
                            dataset=self.name, group=group.name
                        )
                    )
                    asset = self._load_asset(
                        unloaded_asset,
                        key=name,
                        config=data_asset_config,
                        owner=group,
                        putative=putative,
                    )
                    self.groups[group.name].data_assets[name] = asset
            case "sample":
                for sample in self.samples.values():
                    unloaded_asset = Path(
                        str(location_template).format(
                            dataset=self.name, sample=sample.name
                        )
                    )
                    asset = self._load_asset(
                        unloaded_asset,
                        key=name,
                        config=data_asset_config,
                        owner=sample,
                        putative=putative,
                    )
                    self.samples[sample.name].data_assets[name] = asset

    ################################
    # Data Assets Accessors
    ################################

    def get_assets(self, filter_to: list[str] = None) -> list[DataAsset]:
        # extract assets at dataset level
        assets_found: list[DataAsset] = list(self.data_assets.values())

        # extract assets at group level
        for group in self.groups.values():
            assets_found.extend(list(group.data_assets.values()))

        # extract assets at sample level
        for sample in self.samples.values():
            assets_found.extend(list(sample.data_assets.values()))

        if filter_to is not None:
            # TODO: Revisit the ignore here once this issue is resolved. ref: https://github.com/python/mypy/issues/12682
            assets_found = list(filter(lambda x: x.key in filter_to, assets_found))  # type: ignore

        return assets_found

    def compile_multiqc_data(self, data_asset_keys: list[str] = None):
        assets = self.get_assets(filter_to=data_asset_keys)
        return multiqc_run_to_dataframes([asset.path for asset in assets])


def multiqc_run_to_dataframes(paths: list[Path]) -> dict:
    try:
        mqc_ret = multiqc.run(
            analysis_dir=paths,
            quiet=True,
            no_ansi=True,
            no_report=True,
            no_data_dir=True,
            plots_interactive=True,  # ensure data is robustly populated (otherwise flat plots result in missing extractable data)
            # module=[
            #     module.lower() for module in mqc_target["mqc_modules"]
            # ],  # module names here are always lowercase
        )
    except SystemExit:
        raise ValueError(
            f"MultiQC tried to sys.exit. This was given a multiqc.run using these paths: {paths}"
        )

    # extract and set general stats
    general_stats_data = get_general_stats(mqc_ret)
    general_stats = dict()
    for module, data in general_stats_data.items():
        general_stats[module] = pd.DataFrame(
            data
        ).T  # Transpose for consistency with plots dataframes, a samples are the index

    # extract and set plot data
    df_mqc = format_plots_as_dataframe(mqc_ret)

    plots: dict[str, dict[str, pd.DataFrame]] = defaultdict(dict)
    for gb_name, df_gb in df_mqc.groupby(level=[0, 1], axis="columns"):
        # clean dataframe
        # remove row index (redundant with entity ownership)
        # df_gb.reset_index(inplace=True, drop=True)

        # remove top two levels of multindex (these are redundant with gb_name)
        df_gb = df_gb.droplevel(level=[0, 1], axis="columns")

        # Second level of name may indicate subplot (e.g. adapter columns)
        # clean name if the second level is just a copy of the first
        if gb_name[0] == gb_name[1]:
            gb_name = gb_name[0]
        else:
            gb_name = f"{gb_name[0]}:Subplot:{gb_name[1].replace(gb_name[0],'')}"

        # Peel off module name
        module_name, *plot_name = gb_name.split(":")
        # Clean up plot string
        plot_name = ":".join(plot_name).strip()

        plots[module_name][plot_name] = df_gb.dropna(
            how="all"
        )  # remove rows with no data (Usually related to 'samples' not present in all plots, e.g. reads as samples vs samples as samples)

    return {"plots": plots, "general_stats": general_stats}


@dataclass
class Sample:
    name: str
    data_assets: DataAssetDict = field(default_factory=dict, repr=False)

    def compile_multiqc_data(self, data_asset_keys: list[str] = None):
        # TODO: Low priority -> Rework using filter
        if data_asset_keys is None:
            paths = [asset.path for asset in self.data_assets.values()]
        else:
            paths = [
                asset.path
                for asset in self.data_assets.values()
                if asset.key in data_asset_keys
            ]
        return multiqc_run_to_dataframes(paths)


@dataclass
class Group:
    name: str
    data_assets: DataAssetDict = field(default_factory=dict, repr=False)

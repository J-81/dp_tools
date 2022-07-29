from dataclasses import dataclass, field
from pathlib import Path
from string import Formatter
import uuid
from typing import TypedDict
import logging

import pandas as pd

log = logging.getLogger(__name__)


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

    @staticmethod
    def _create_asset(
        asset: Path, key: str, config: dict, owner: ExperimentalEntity
    ) -> DataAsset:
        if "*" in asset.name:
            [asset] = asset.parent.glob(asset.name)
        assert asset.exists(), f"Failed to load asset at path '{asset}'"
        return DataAsset(key=key, path=asset, config=config, owner=owner)

    # TODO: dict -> better typehint via typeddict
    def load_data_asset(self, data_asset_config: dict, root_dir: Path, name: str):
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
                asset = self._create_asset(
                    unloaded_asset, key=name, config=data_asset_config, owner=self
                )
                self.data_assets[name] = asset
            case "group":
                for group in self.groups:
                    unloaded_asset = Path(
                        str(location_template).format(dataset=self.name, group=group)
                    )
                    asset = self._create_asset(
                        unloaded_asset, key=name, config=data_asset_config, owner=group
                    )
                    self.groups[group].data_assets[name] = asset
            case "sample":
                for sample in self.samples:
                    unloaded_asset = Path(
                        str(location_template).format(dataset=self.name, sample=sample)
                    )
                    asset = self._create_asset(
                        unloaded_asset, key=name, config=data_asset_config, owner=sample
                    )
                    self.samples[sample].data_assets[name] = asset

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


@dataclass
class Sample:
    name: str
    data_assets: DataAssetDict = field(default_factory=dict, repr=False)


@dataclass
class Group:
    name: str
    data_assets: DataAssetDict = field(default_factory=dict, repr=False)

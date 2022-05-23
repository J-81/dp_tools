from ast import Assert
from dataclasses import dataclass, field
import abc
import datetime
import enum
import hashlib
from pathlib import Path
import uuid
from typing import (
    Dict,
    List,
    Literal,
    Optional,
    OrderedDict,
    Protocol,
    Set,
    Tuple,
    TypedDict,
    Union,
    runtime_checkable,
)
import logging
from click import File

import pandas as pd
import multiqc

from dp_tools.core.model_commons import strict_type_checks
from dp_tools.core.utilites.multiqc_tools import (
    format_plots_as_dataframe,
    get_general_stats,
)

log = logging.getLogger(__name__)


def get_id():
    uuid4 = uuid.uuid4()
    return str(uuid4)


#########################################################################
# DATAFILE
#########################################################################
@dataclass
class DataFile:
    """A class for keeping track of files"""

    path: Path
    metadata: dict = field(
        default_factory=dict
    )  # used to incorporate to additional data asset metadata like publish status etc.
    md5sum: str = field(init=False, repr=False)  # compute on ingestion
    # id: str = field(default_factory=get_id)
    dummy_md5sum: bool = field(
        default=False
    )  # used to replace md5sum of contents with md5sum of the file name
    check_exists: bool = field(default=True)

    def __post_init__(self):
        log.debug(f"Initiating DataFile with path: {self.path}")
        if self.check_exists:
            if not self.path.is_file():
                raise FileNotFoundError(
                    f"Cannot create DataFile with non-existent file: {self.path}"
                )
        if self.dummy_md5sum:
            # generate md5sum based on the path/file name only
            # this computes extremely fast and is advised for testing
            self.md5sum = f"DUMMY:{self._compute_md5sum(self.path.name.encode())}"
            self._md5sum_cached = True
        else:
            self._md5sum_cached = False

        # finally enforce types check
        strict_type_checks(self, exceptions=["md5sum"])

    def compute_md5sum(self):
        if not self._md5sum_cached:
            self.md5sum = self._compute_md5sum(self.path.open("rb").read())
        return self.md5sum

    def _compute_md5sum(self, contents):
        return hashlib.md5(contents).hexdigest()

    def __str__(self):
        return f"{self.__class__.__name__}(path={self.path})"


@dataclass
class DataDir:
    """A class for keeping track of directories
    Both multiQC and RSEM include examples of things better tracked as directories; however, more granual loading can be achieved with DataFile.
    """

    path: Path
    metadata: dict = field(
        default_factory=dict
    )  # used to incorporate to additional data asset metadata like publish status etc.
    check_exists: bool = field(default=True)

    def __post_init__(self):
        log.debug(f"Initiating DataDir with path: {self.path}")
        if self.check_exists:
            if not self.path.is_dir():
                raise FileNotFoundError(
                    f"Cannot create DataDir with non-existent directory: {self.path}"
                )

        # finally enforce types check
        strict_type_checks(self)

    def __str__(self):
        return f"{self.__class__.__name__}(path={self.path})"


#########################################################################
# FUNCTIONAL MIXINS
# THESE MAY IMPART SHARED FUNCTIONALITY ACROSS LEVES
# E.G. Attaching components should work on datasets and samples
#########################################################################
# a mixin class
class CanAttachComponents:
    def _is_component(self, putative_component):
        """Check is the component is a component
        Specifically, check if the base is a BaseComponent
        """
        return isinstance(getattr(putative_component, "base", None), BaseComponent)

    @property
    def all_components(self) -> Dict[str, object]:
        """Return a dictionary of components, including any empty slots"""
        components = dict()
        for attr, value in self.__dict__.items():
            # check the base, all component bases should be a BaseComponent
            # if 'base' exists, check if it is a component
            if self._is_component(value):
                components[attr] = value
        return components

    @property
    def _str_all_components(self) -> str:
        """Return a dictionary of components, including any empty slots"""
        s = (
            "Components:\n\t"
            + " \n\t".join(
                [
                    f"Attribute '{attr}' --> '{comp.__class__.__name__}'"
                    for attr, comp in self.all_components.items()
                ]
            )
            if self.all_components
            else "No Components. (including empty component slots)"
        )
        return s

    @property
    def components(self) -> List[object]:
        """Return a list of components, EXCLUDING any empty slots"""
        components = list()
        for attr, value in self.__dict__.items():
            # check the base, all component bases should be a BaseComponent
            # if 'base' exists, check if it is a component
            if self._is_component(value) and not isinstance(value, EmptyComponent):
                components.append(value)
        return components

    def attach_component(self, component, attr):
        # ensure this is an expected component by name
        # TODO: replace this with only checking component attrs,
        # perhaps by checking if default is a component since EmptyComponent will return True
        assert hasattr(
            self, attr
        ), "Sample does not specify component with slot: {attr}"

        # ensure this is a component
        if not self._is_component(component):
            raise TypeError(
                "Failure: Cannot attach unless object 'base' is a BaseComponent"
            )

        # warn if component already attached and not an EmptyComponent
        if not isinstance(getattr(self, attr), EmptyComponent):
            log.warning(f"Overwriting existing component in slot: {attr}")

        # attach component
        log.debug(f"attaching {attr} to {self.name}")
        setattr(self, attr, component)
        log.debug(f"post attaching {attr} to {self.name}")

        # associate component with sample
        log.debug(f"attaching {self.name} to {attr}")
        component._attach_entity(entity=self, entity_attr=attr)
        log.debug(f"post attaching {self.name} to {attr}")


# Dict[attr: str, 'entity': Union['TemplateDataset', 'TemplateSample' ]
class AttachEntity(TypedDict):
    attr: str
    entity: Union["TemplateDataset", "TemplateSample"]


# mixin for attaching samples or datasets or datasystems
class CanAttachEntity:
    entities: Dict[str, AttachEntity] = None  # type: ignore
    entityType: str = None  # type: ignore

    def _attach_entity(
        self, entity: Union["TemplateDataset", "TemplateSample"], entity_attr: str
    ):
        # create dict if first time for this component
        if not self.entities:
            self.entities = dict()
            if isinstance(entity, TemplateSample):
                self.entityType = "sample"
            elif isinstance(entity, TemplateDataset):
                self.entityType = "dataset"
            self._entityClass = entity.__class__

        # block mixing entity classes
        if not isinstance(entity, self._entityClass):
            raise ValueError(
                f"Cannot attach entity type: {entity.__class__} after already attaching to entity class: {self._entityClass}"
            )
        self.entities[entity.name] = {"entity": entity, "attr": entity_attr}

    @property
    def dataset(self) -> "TemplateDataset":
        match self.entityType:
            case "sample":
                dataset = list(self.entities.values())[0]["entity"].dataset
            case "dataset":
                dataset = list(self.entities.values())[0]["entity"]
        return dataset


#########################################################################
# COMPONENTS
#########################################################################
# type hint definitions
SampleLevelDict = OrderedDict[  # samples level
    "str",  # denotes sample name
    OrderedDict["str", float],  # key:value level, e.g. 'percent_gc': 52.0
]

GeneralStatsDict = Dict["str", SampleLevelDict]  # modules level  # denotes module name

GeneralStatsDictForComponent = OrderedDict[  # samples level
    "str",  # denotes sample name
    OrderedDict["str", float],  # key:value level, e.g. 'percent_gc': 52.0
]

DataAsset = Union[DataDir, DataFile]

# model dict types
TopLevelMQC = Dict[str, "ModuleLevelMQC"]


class ModuleLevelMQC(TypedDict):
    General_Stats: Dict[str, Union[int, float]]  # holds general stat values
    Plots: Dict[str, pd.DataFrame]  # holds plotname: dataframe


# define type hint
class MQCTargetDict(TypedDict):
    mqc_modules: List[str]  # Module names to run
    target: DataAsset


@dataclass
class BaseComponent:
    """Class for keeping track of abstract components like Reads, Alignments, and Counts"""

    # id: str = field(default_factory=get_id)
    description: str = field(default="Null")
    created: datetime.datetime = field(default_factory=datetime.datetime.now)


class TemplateComponent(abc.ABC, CanAttachEntity):
    # TODO: refine this type hint
    _mqc_data: dict
    _mqcData: Union[Dict, bool]  # only set to false for EmptyComponent
    # _mqc_modules_limit: List[str] = [
    #     "FastQC"
    # ]  # optional but speeds up multiqc processing time slightly and can constrain modules if desired

    def __post_init__(self):
        log.debug(f"Initiating Component with class: {self.__class__.__name__}")
        strict_type_checks(self)

    # TODO: refactor, attached entity related portion should be owned by mixin
    def __repr__(self):
        attachments_str = (
            "Attached to Entity:\n\t"
            + "\n\t".join(
                [
                    f'{entity_name} <--AS-- {entity_dict["attr"]}'
                    for entity_name, entity_dict in self.entities.items()
                ]
            )
            if self.entities
            else ""
        )
        data_str = (
            "Has Data Assets (listed by attribute name):\n\t"
            + "\n\t".join(
                [
                    f"'{attr}' :{data}"
                    for attr, data in self.__dict__.items()
                    if any([isinstance(data, DataDir), isinstance(data, DataFile)])
                ]
            )
            if self.__dict__
            else ""
        )
        return f"""
Component: class:{self.__class__.__name__}
{data_str}
{attachments_str}
        """

    def __str__(self):
        return f"""{self.__class__.__name__}"""

    @property
    def all_data_assets(self):
        assets = list()
        for value in self.__dict__.values():
            if any([isinstance(value, DataFile), isinstance(value, DataDir)]):
                assets.append(value)
        return assets

    @property
    def _mqc_targets(self) -> List[MQCTargetDict]:
        """A list of all components that have a multiQC module available for parsing
        This is determined using the metadata field, specifically 'mqc_parse'
        """
        all_targets = list()
        # iterate through all dataclass fields
        for field_attr, field in self.__dataclass_fields__.items():
            # check that data asset is field-metadata marked
            if mqc_modules := field.metadata.get("mqc_parse"):
                target: MQCTargetDict = {
                    "mqc_modules": mqc_modules,
                    "target": getattr(self, field_attr),
                }
                all_targets.append(target)
        return all_targets

    @property
    def mqcData(self) -> dict[str, ModuleLevelMQC]:
        if getattr(self, "_mqcData", None) == None:
            # extract mqcData
            self._mqcData: dict[str, ModuleLevelMQC] = dict()
            for mqc_target in self._mqc_targets:
                log.debug(f"Extracting multiqc data for {mqc_target}")
                self.__extract_mqcData(mqc_target)
        return self._mqcData

    def __extract_mqcData(self, mqc_target: MQCTargetDict):
        """Populates data for a specific target for
        Includes both general stats and plots.
        Can include multiple modules as designated in the mqc_target dictionary.
        """
        assert len(
            mqc_target["mqc_modules"]
        ), "Only one module per data asset currently supported"
        # set module_name variable for plot extraction usage
        module_name = mqc_target["mqc_modules"][0]
        # run MQC
        # skip if the target doesn't exist (e.g. an unattached data asset)
        if not mqc_target["target"]:
            log.warning(
                f"{mqc_target['target']} designated as a mqc target, but data asset was missing"
            )
            return
        mqc_ret = multiqc.run(
            analysis_dir=[mqc_target["target"].path],
            no_report=True,
            no_data_dir=True,
            plots_interactive=True,  # ensure data is robustly populated (otherwise flat plots result in missing extractable data)
            module=[
                module.lower() for module in mqc_target["mqc_modules"]
            ],  # module names here are always lowercase
        )

        # extract and set general stats
        general_stats = get_general_stats(mqc_ret)
        if not self._mqcData.get(module_name):
            self._mqcData[module_name]: ModuleLevelMQC = dict()  # type: ignore
        # handle modules with no general stats
        if general_stats == {}:
            log.info(
                f"MQC extraction: No general stats in found after running: {mqc_target['mqc_modules']}"
            )
            self._mqcData[module_name]["General_Stats"] = False
        else:
            for module_name, samplewise_dict in general_stats.items():
                # init dict if needed
                # TODO: replace with default dict on toplevel init
                if not self._mqcData.get(module_name):
                    self._mqcData[module_name]: ModuleLevelMQC = dict()  # type: ignore

                # assert only one mqc-sample generated
                assert (
                    len(list(samplewise_dict.values())) == 1
                ), "For unshared logs, this is a true issue.  This will break on shared logs"
                # remove sample layer on assignment (instead sample is implicit in component entity)
                self._mqcData[module_name]["General_Stats"] = list(
                    samplewise_dict.values()
                )[0]

        # extract and set plot data
        df_mqc = format_plots_as_dataframe(mqc_ret)

        for gb_name, df_gb in df_mqc.groupby(level=[0, 1], axis="columns"):
            # clean dataframe
            # remove row index (redundant with entity ownership)
            df_gb.reset_index(inplace=True, drop=True)

            # remove top two levels of multindex (these are redundant with gb_name)
            df_gb = df_gb.droplevel(level=[0, 1], axis="columns")

            # init dict if needed
            # TODO: replace with default dict on toplevel init
            if not self._mqcData.get(module_name):
                self._mqcData[module_name]: ModuleLevelMQC = dict()  # type: ignore
            if not self._mqcData[module_name].get("Plots"):
                self._mqcData[module_name]["Plots"] = dict()

            # clean name
            # first remove module name as that is a higher key making it redundant
            remove_substr = f"{module_name}: "
            gb_name = (
                gb_name[0].replace(remove_substr, ""),
                gb_name[1].replace(remove_substr, ""),
            )

            # Second level of name may indicate subplot (e.g. adapter columns)
            # clean name if the second level is just a copy of the first
            if gb_name[0] == gb_name[1]:
                gb_name = gb_name[0]
            else:
                gb_name = f"{gb_name[0]}:Subplot:{gb_name[1].replace(gb_name[0],'')}"

            # remove sample layer on assignment (instead sample is implicit in component entity)
            self._mqcData[module_name]["Plots"][gb_name] = df_gb


@dataclass(eq=False)
class EmptyComponent(TemplateComponent):
    """Class representing an empty component"""

    base: BaseComponent = BaseComponent(description="This slot is empty")
    _mqcData: bool = False  # only set to false for EmptyComponent


#########################################################################
# DATASYSTEM
#########################################################################
@dataclass
class BaseDataSystem:
    """Abstract class for a dataset
    NOTE: Within this system:
        a dataset is associated with a specific assay
        a dataSYSTEM is the container for multiple datasets
    """

    name: str
    datasets: Dict[str, "BaseDataset"] = field(default_factory=dict, repr=False)


class TemplateDataSystem(abc.ABC):
    """This abstract base class should serve as a template for new data systems"""

    allowed_dataset_classes: list

    @property
    def name(self):
        """Alias name to from compositioned base to implementation"""
        return self.base.name

    @property
    def dataset(self):
        """Convenience function when only one dataset is attached"""
        if len(self.base.datasets) == 0:
            return None
        elif len(self.base.datasets) == 1:
            return list(self.base.datasets.values())[0]
        else:
            raise ValueError(
                f"This datasystem has multiple datasets. Use 'datasets' or 'all_datasets' instead to indicate specific dataset"
            )

    @property
    def datasets(self):
        return self.base.datasets

    @property
    def all_datasets(self):
        """returns a set of all datasets"""
        return set(self.datasets.values())

    @property
    def all_samples(self):
        """returns a set of all samples"""
        all_samples = set()
        for dataset in self.all_datasets:
            all_samples = all_samples.union(set(dataset.samples.values()))
        return all_samples

    @property
    def all_components(self):
        """returns a set of all components"""
        all_components = set()
        # get dataset wise components
        for dataset in self.all_datasets:
            all_components = all_components.union(dataset.all_components.values())

        # get sample wise components
        for sample in self.all_samples:
            all_components = all_components.union(sample.all_components.values())

        return all_components

    def attach_dataset(self, dataset):
        if dataset.name in self.datasets.keys():
            log.warning(f"Overwriting pre-existing dataset: {dataset.name}")
        self.datasets[dataset.name] = dataset

        # associated dataset to datasystem
        dataset.dataSystem = self

    def __str__(self):
        dataset_sub_s = "\n\t".join(
            [dataset.name for dataset in self.datasets.values()]
        )
        return f"""
DataSystem (class:{self.__class__.__name__}): '{self.name}'
Has Datasets: 
    {dataset_sub_s}
        """


#########################################################################
# DATASET
#########################################################################
@dataclass
class BaseDataset:
    """Abstract class for a dataset
    NOTE: Within this system:
        a dataset is associated with a specific assay
        a dataSYSTEM is the container for multiple datasets
    """

    # id: str = field(default_factory=get_id)
    name: str = field(default="orphan")
    samples: Dict[str, "BaseSample"] = field(default_factory=dict, repr=False)
    config: dict = field(default_factory=dict, repr=False)


# for subclassing
class TemplateDataset(abc.ABC, CanAttachComponents):
    dataSystem: TemplateDataSystem = None

    @property
    @abc.abstractproperty
    def expected_sample_class(self):
        return self.expected_sample_class

    @property
    def samples(self):
        return self.base.samples

    @property
    def name(self):
        return self.base.name

    @property
    def config(self):
        return self.base.config

    # TODO: add type, how to type hint a class in general
    def attach_sample(self, sample):
        if not isinstance(sample, self.expected_sample_class):
            raise TypeError(
                f"Improper sample attachment class: expected-{self.expected_sample_class} got-{type(sample)}"
            )
        if sample.name in self.samples.keys():
            log.warning(f"Overwriting pre-existing sample: {sample.name}")
        # attach sample to dataset
        self.samples[sample.name] = sample
        # attach dataset to sample
        sample.base.dataset = self

    # extends mixin CanAttachComponents
    @property
    def all_non_empty_components_recursive(self):
        """returns a set of all components"""
        all_components = set()
        # get dataset attached components
        all_components = all_components.union(
            [
                comp
                for comp in self.all_components.values()
                if not isinstance(comp, EmptyComponent)
            ]
        )

        # get sample wise components
        for sample in self.samples.values():
            all_components = all_components.union(
                [
                    comp
                    for comp in sample.all_components.values()
                    if not isinstance(comp, EmptyComponent)
                ]
            )

        return all_components

    @property
    def all_data_assets(self):
        """returns a set of all components"""
        assets = list()
        for component in self.all_non_empty_components_recursive:
            assets.extend(component.all_data_assets)
        return assets

    def __str__(self):
        return f"""
Dataset (class:{self.__class__.__name__}): '{self.name}'
Has {len(self.samples)} Samples (class:{self.expected_sample_class.__name__})
{self._str_all_components}
Is Part of DataSystem: 
    '{self.dataSystem.name}'
        """

    # TODO: refactor shared parts of getMQCPlots and getMQCDataFrame
    def getMQCPlots(self, sample_component: str, mqc_module: str) -> list[str]:
        # check if sample requested component exists for all samples
        has_component = {
            sample.name: getattr(sample, sample_component, False)
            for sample in self.samples.values()
        }
        assert all(
            has_component.values()
        ), f"At least one sample does not have the requested component ({sample_component})!: {has_component}"

        # check if all samples have the requested mqc module
        has_module = {
            sample.name: getattr(sample, sample_component).mqcData.get(mqc_module)
            for sample in self.samples.values()
        }
        assert all(
            has_module.values()
        ), f"At least one sample does not have the requested multiQC module ({mqc_module})!: {has_module}"

        # return list of plots
        plots_per_sample: List[Set] = [
            set(
                getattr(sample, sample_component)
                .mqcData.get(mqc_module)["Plots"]
                .keys()
            )
            for sample in self.samples.values()
        ]
        # find the intersection of sample wise sets of plots
        plots_in_all_samples = (
            set.intersection(*plots_per_sample) if plots_per_sample else set()
        )  # guard against empty list
        plots_in_all_samples.add("general_stats")
        return list(plots_in_all_samples)

    def getMQCDataFrame(
        self, sample_component: str, mqc_module: str, mqc_plot: str = None
    ) -> pd.DataFrame:
        """Generates a single dataframe composed of all samples for the requested component, mqc module, and plot.
        If a mqc_plot is not specificed, a list of available plots is listed.

        :param sample_component: The sample component attribute name
        :type sample_component: str
        :param mqc_module: The target components multiQC module name
        :type mqc_module: str
        :param mqc_plot: The target components plot as generated by the requested multiQC module
        :type mqc_plot: str
        :return: A dataframe with the index composed of sample names and the columns consistent with the original sample-wise multiQC plot dataframe
        :rtype: pd.DataFrame
        """
        # check if sample requested component exists for all samples
        has_component = {
            sample.name: getattr(sample, sample_component, False)
            for sample in self.samples.values()
        }
        assert all(
            has_component.values()
        ), f"At least one sample does not have the requested component ({sample_component})!: {has_component}"

        # check if all samples have the requested mqc module
        has_module = {
            sample.name: getattr(sample, sample_component).mqcData.get(mqc_module)
            for sample in self.samples.values()
        }
        assert all(
            has_module.values()
        ), f"At least one sample does not have the requested multiQC module ({mqc_module})!: {has_module}"

        # return table of general stats instead
        if mqc_plot == "general_stats":
            dict_format = {
                s.name: getattr(s, sample_component).mqcData[mqc_module][
                    "General_Stats"
                ]
                for s in self.samples.values()
            }
            return pd.DataFrame(dict_format).T

        # check if all samples have the requested mqc plot
        # need to explicitly check if the returned is a dataframes
        # as dataframe don't eval to True by default (ambiguous)
        # additionally checking if it is None does a cell by cell check and returns a Dataframe too
        has_plot = {
            sample.name: isinstance(
                getattr(sample, sample_component)
                .mqcData.get(mqc_module)["Plots"]
                .get(mqc_plot),
                pd.DataFrame,
            )
            for sample in self.samples.values()
        }
        assert all(
            has_plot.values()
        ), f"At least one sample does not have the requested plot ({mqc_plot}) in module ({mqc_module})!: {has_plot}"

        return pd.concat(
            [
                getattr(s, sample_component)
                .mqcData[mqc_module]["Plots"][mqc_plot]
                .set_index(pd.Index([s.name]))
                for s in self.samples.values()
            ]
        )


#########################################################################
# SAMPLES
#########################################################################
@dataclass
class BaseSample:
    """Abstract class for samples"""

    name: str
    # id: str = field(default_factory=get_id)
    dataset: Union[None, BaseDataset] = field(default=None)


class TemplateSample(abc.ABC, CanAttachComponents):
    # used properties to alias base attributes as needed
    @property
    def name(self):
        return self.base.name

    @property
    def dataset(self):
        return self.base.dataset

    def __str__(self):
        return f"""
Sample (class:{self.__class__.__name__}): '{self.name}'
{self._str_all_components}
Is Part of Dataset: 
    '{self.dataset.name if self.dataset else "None"}'
        """


############################################################################################
# GLDS SPECIFIC
############################################################################################
@dataclass(eq=False)
class GLDSDataSystem(TemplateDataSystem):

    base: BaseDataSystem = field(repr=False)

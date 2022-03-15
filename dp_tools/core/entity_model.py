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
    Protocol,
    Tuple,
    TypedDict,
    Union,
    runtime_checkable,
)
import logging

import pandas as pd

from dp_tools.core.model_commons import strict_type_checks

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)


def get_id():
    uuid4 = uuid.uuid4()
    return str(uuid4)


#########################################################################
# DATAFILE
#########################################################################
@dataclass
class DataFile:
    """ A class for keeping track of files """

    path: Path
    md5sum: str = field(init=False)  # compute on ingestion
    # id: str = field(default_factory=get_id)
    dummy_md5sum: bool = field(
        default=False
    )  # used to replace md5sum of contents with md5sum of the file name

    def __post_init__(self):
        assert self.path.is_file()
        if self.dummy_md5sum:
            # generate md5sum based on the path/file name only
            # this computes extremely fast and is advised for testing
            self.md5sum = f"DUMMY:{self._compute_md5sum(self.path.name.encode())}"
        else:
            self.md5sum = self._compute_md5sum(self.path.open("rb").read())

        # finally enforce types check
        strict_type_checks(self)

    def _compute_md5sum(self, contents):
        return hashlib.md5(contents).hexdigest()


@dataclass
class DataDir:
    """ A class for keeping track of directories
    Both multiQC and RSEM include examples of things better tracked as directories; however, more granual loading can be achieved with DataFile.
    """

    path: Path

    def __post_init__(self):
        assert self.path.is_dir()
        # finally enforce types check
        strict_type_checks(self)


#########################################################################
# FUNCTIONAL MIXINS
# THESE MAY IMPART SHARED FUNCTIONALITY ACROSS LEVES
# E.G. Attaching components should work on datasets and samples
#########################################################################
# a mixin class
class CanAttachComponents:
    def _is_component(self, putative_component):
        """ Check is the component is a component
        Specifically, check if the base is a BaseComponent
        """
        return isinstance(getattr(putative_component, "base", None), BaseComponent)

    @property
    def all_components(self) -> Dict[str, object]:
        """ Return a dictionary of components, including any empty slots """
        components = dict()
        for attr, value in self.__dict__.items():
            # check the base, all component bases should be a BaseComponent
            # if 'base' exists, check if it is a component
            if self._is_component(value):
                components[attr] = value
        return components

    @property
    def components(self) -> List[object]:
        """ Return a list of components, EXCLUDING any empty slots """
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
        log.debug(f"attaching {component} to {self}")
        setattr(self, attr, component)
        log.debug(f"post attaching {component} to {self}")

        # associate component with sample
        log.debug(f"attaching {self.name} to {component}")
        component._attach_entity(entity=self, entity_attr=attr)
        log.debug(f"post attaching {self.name} to {component}")

# Dict[attr: str, 'entity': Union['TemplateDataset', 'TemplateSample' ]
class AttachEntity(TypedDict):
    attr: str
    entity: Union['TemplateDataset', 'TemplateSample']

# mixin for attaching samples or datasets or datasystems
class CanAttachEntity:
    entities: Dict[str, AttachEntity] = None
    entityType: str = None

    def _attach_entity(self, entity: Union['TemplateDataset', 'TemplateSample' ], entity_attr: str):
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
            raise ValueError(f"Cannot attach entity type: {entity.__class__} after already attaching to entity class: {self._entityClass}")
        self.entities[entity.name] = {'entity':entity, 'attr': entity_attr}

    @property
    def dataset(self) -> 'TemplateDataset':
        match self.entityType:
            case "sample":
                dataset = list(self.entities.values())[0]["entity"].dataset
            case "dataset":
                dataset = list(self.entities.values())[0]['entity']
        return dataset


#########################################################################
# COMPONENTS
#########################################################################
@dataclass
class BaseComponent:
    """ Class for keeping track of abstract components like Reads, Alignments, and Counts """

    # id: str = field(default_factory=get_id)
    description: str = field(default="Null")
    created: datetime.datetime = field(default_factory=datetime.datetime.now)


class TemplateComponent(abc.ABC):
    def __post_init__(self):
        strict_type_checks(self)


@dataclass(eq=False)
class EmptyComponent(TemplateComponent):
    """ Class representing an empty component """

    base: BaseComponent = BaseComponent(description="This slot is empty")


#########################################################################
# DATASYSTEM
#########################################################################
@dataclass
class BaseDataSystem:
    """ Abstract class for a dataset
    NOTE: Within this system: 
        a dataset is associated with a specific assay
        a dataSYSTEM is the container for multiple datasets
    """

    name: str
    datasets: Dict[str, "BaseDataset"] = field(default_factory=dict, repr=False)


class TemplateDataSystem(abc.ABC):
    """ This abstract base class should serve as a template for new data systems """

    allowed_dataset_classes: list

    @property
    def name(self):
        """ Alias name to from compositioned base to implementation """
        return self.base.name

    @property
    def dataset(self):
        """ Convenience function when only one dataset is attached """
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
        """ returns a set of all datasets """
        return set(self.datasets.values())

    @property
    def all_samples(self):
        """ returns a set of all samples """
        all_samples = set()
        for dataset in self.all_datasets:
            all_samples = all_samples.union(set(dataset.samples.values()))
        return all_samples

    @property
    def all_components(self):
        """ returns a set of all components """
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
        return f"""
        DataSystem ({self.__class__}): {self.name}
        Has Datasets: {self.datasets}
        """

#########################################################################
# DATASET
#########################################################################
@dataclass
class BaseDataset:
    """ Abstract class for a dataset
    NOTE: Within this system: 
        a dataset is associated with a specific assay
        a dataSYSTEM is the container for multiple datasets
    """

    # id: str = field(default_factory=get_id)
    name: str = field(default="orphan")
    samples: Dict[str, "BaseSample"] = field(default_factory=dict, repr=False)


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
        """ returns a set of all components """
        all_components = set()
        # get dataset attached components
        all_components = all_components.union([comp for comp in self.all_components.values() if not isinstance(comp, EmptyComponent)])

        # get sample wise components
        for sample in self.samples.values():
            all_components = all_components.union([comp for comp in sample.all_components.values() if not isinstance(comp, EmptyComponent)])

        return all_components

    def __str__(self):
        return f"""
        Dataset ({self.__class__}): {self.name}
        Has {len(self.samples)} Samples ({self.expected_sample_class})
        Dataset Components: {self.all_components}
        Is Part of DataSystem: {self.dataSystem.name}
        """


#########################################################################
# SAMPLES
#########################################################################
@dataclass
class BaseSample:
    """ Abstract class for samples """

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
        Sample ({self.__class__}): {self.name}
        Sample Components: {self.all_components}
        Is Part of Dataset: {self.dataset}
        """

############################################################################################
# GLDS SPECIFIC
############################################################################################
@dataclass(eq=False)
class GLDSDataSystem(TemplateDataSystem):

    base: BaseDataSystem = field(repr=False)

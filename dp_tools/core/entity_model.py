from ast import Assert
from dataclasses import dataclass, field
import abc
import datetime
import enum
import hashlib
from pathlib import Path
import uuid
from typing import Dict, List, Literal, Optional, Protocol, Tuple, Union, runtime_checkable
import logging
from dp_tools.core.check_model import Flag, FlagCode

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


class MustValidate:
    """ A mixin to add in a common validation method """

    @abc.abstractmethod
    def validate(self) -> List[Flag]:
        ...

    def _strict_validate(self):
        flags = self.validate()
        # check that all returned objects are indeed Flag objects
        assert all([isinstance(flag, Flag) for flag in flags])
        # log a warning if no flags are returned
        if flags == []:
            log.warning(f"This object generated no validation flags: {self}")
        return flags


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
        component._attach_entity(self)
        log.debug(f"post attaching {self.name} to {component}")


# mixin for attaching samples or datasets or datasystems
class CanAttachEntity:
    entities: dict = dict()

    def _attach_entity(self, entity):
        self.entities[entity.name] = entity


#########################################################################
# COMPONENTS
#########################################################################
@dataclass
class BaseComponent:
    """ Class for keeping track of abstract components like Reads, Alignments, and Counts """

    # id: str = field(default_factory=get_id)
    description: str = field(default="Null")
    created: datetime.datetime = field(default_factory=datetime.datetime.now)


class TemplateComponent(abc.ABC, MustValidate):
    def __post_init__(self):
        strict_type_checks(self)


@dataclass(eq=False)
class EmptyComponent(TemplateComponent):
    """ Class representing an empty component """

    base: BaseComponent = BaseComponent(description="This slot is empty")

    def validate(self) -> List[Tuple]:
        flag = Flag(self, description="This is an empty component.", code=FlagCode.INFO)
        return [flag]


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

    def validate(self):
        log.info("Running validation at data system level")
        dataset_flags = {
            dataset: dataset._strict_validate() for dataset in self.all_datasets
        }  # return of validate should be flags
        sample_flags = {
            sample: sample._strict_validate() for sample in self.all_samples
        }
        # return flags solely from non-Empty components
        component_flags = {
            component: component._strict_validate()
            for component in self.all_components
            if not isinstance(component, EmptyComponent)
        }
        if not dataset_flags:
            log.warning(f"No datasets found nor validated in entire dataSystem")
        if not sample_flags:
            log.warning(f"No samples found nor validated in entire dataSystem")
        if not component_flags:
            log.warning(f"No components found nor validated in entire dataSystem")

        # TODO: iterate through flags and check if anything should raise an exitcode
        # maybe hand these flags off to a halt assessor?
        log.debug(
            f"Handing these flags off to TBD: dataset: {dataset_flags}, sample: {sample_flags}, component_flags: {component_flags}"
        )

        self.validation_report = dict()
        self.validation_report["dataset"] = dataset_flags
        self.validation_report["sample"] = sample_flags
        self.validation_report["component"] = component_flags

        # validation report should not perform actions itself
        # actions might include report generation and raising an Exception


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
class TemplateDataset(abc.ABC, MustValidate):
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


#########################################################################
# SAMPLES
#########################################################################
@dataclass
class BaseSample:
    """ Abstract class for samples """

    name: str
    # id: str = field(default_factory=get_id)
    dataset: Union[None, BaseDataset] = field(default=None)


class TemplateSample(abc.ABC, MustValidate):
    # used properties to alias base attributes as needed
    @property
    def name(self):
        return self.base.name

    @property
    def dataset(self):
        return self.base.dataset

############################################################################################
# GLDS SPECIFIC
############################################################################################
@dataclass(eq=False)
class GLDSDataSystem(TemplateDataSystem):

    base: BaseDataSystem = field(repr=False)


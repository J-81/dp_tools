from ast import Assert
from dataclasses import dataclass, field
import abc
import datetime
import enum
import hashlib
from pathlib import Path
import uuid
from typing import Dict, List, Literal, Union
import logging

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)


def get_id():
    uuid4 = uuid.uuid4()
    return str(uuid4)


def strict_type_checks(
    obj: object, exceptions: List[str] = None, except_nones: List[str] = None
):
    # set default empty lists
    exceptions = exceptions if exceptions else list()
    except_nones = except_nones if except_nones else list()

    for (name, field_type) in obj.__annotations__.items():
        pass_conditions = list()
        if name in exceptions:
            log.debug(f"Excluding type checking for {name}")
            continue
        if name in except_nones:
            log.debug(f"Allowing 'None' as valid for {name}")
            pass_conditions.append(not obj.__dict__[name])
            continue
        # base type check
        pass_conditions.append(isinstance(obj.__dict__[name], field_type))
        if not any(pass_conditions):
            current_type = type(obj.__dict__[name])
            raise TypeError(
                f"The field `{name}` was assigned by `{current_type}` instead of `{field_type}`"
            )

    log.debug("Strict Type Check is passed successfully")


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

    def attach_component(self, component, attr):
        # ensure this is an expected component by name
        # TODO: replace this with only checking component attrs, 
        # perhaps by checking if default is a component since EmptyComponent will return True
        assert hasattr(
            self, attr
        ), "Sample does not specify component with slot: {attr}"

        # ensure this is a component
        if not self._is_component(component):
            raise TypeError("Failure: Cannot attach unless object 'base' is a BaseComponent")

        # warn if component already attached
        if getattr(self, attr):
            log.warning(f"Overwriting existing component in slot: {attr}")

        # attach component
        log.debug(f"attaching {component} to {self}")
        setattr(self, attr, component)
        log.debug(f"post attaching {component} to {self}")

        # associate component with sample
        log.debug(f"attaching {self.name} to {component}")
        component.attach_entity(self)
        log.debug(f"post attaching {self.name} to {component}")


# mixin for attaching samples or datasets or datasystems
@dataclass
class CanAttachEntity:
    entities: dict = field(default_factory=dict)

    def attach_entity(self, entity):
        self.entities[entity.name] = entity


#########################################################################
# COMPONENTS
#########################################################################
@dataclass
class BaseComponent(CanAttachEntity):
    """ Class for keeping track of abstract components like Reads, Alignments, and Counts """

    # id: str = field(default_factory=get_id)
    description: str = field(default="Null")
    created: datetime.datetime = field(default_factory=datetime.datetime.now)


@dataclass
class EmptyComponent:
    """ Class representing an empty component """

    base: BaseComponent = BaseComponent(description="This slot is empty")


@dataclass
class ReadsComponent:

    base: BaseComponent
    fastqGZ: DataFile
    # id: str = field(default_factory=get_id)
    multiQCZip: Union[DataFile, None] = field(default=None)
    fastqcReportHTML: Union[DataFile, None] = field(default=None)
    fastqcReportZIP: Union[DataFile, None] = field(default=None)
    trimmingReportTXT: Union[DataFile, None] = field(default=None)

    def __post_init__(self):
        strict_type_checks(self)

    @property
    def entities(self):
        return self.base.entities

    def attach_entity(self, *args, **kwargs):
        return self.base.attach_entity(*args, **kwargs)


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


class TemplateDataSystem:
    """ This abstract base class should serve as a template for new data systems """

    allowed_dataset_classes: list

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
            all_samples = all_samples.union(dataset.all_components)
        return all_samples

    @property
    def all_components(self):
        """ returns a set of all components """
        all_components = set()
        # get dataset wise components
        for dataset in self.all_datasets:
            all_components = all_components.union(dataset.all_components)

        # get sample wise components
        for sample in self.all_samples:
            all_components = all_components.union(sample.all_components)

        return all_components

    def attach_dataset(self, dataset):
        if dataset.name in self.datasets.keys():
            log.warning(f"Overwriting pre-existing dataset: {dataset.name}")
        self.datasets[dataset.name] = dataset

    def validate(self):
        log.info("Running validation at data system level")
        dataset_flags = {
            dataset.name: dataset.validate() for dataset in self.all_datasets
        }  # return of validate should be flags
        sample_flags = {sample.name: sample.validate() for sample in self.all_samples}
        component_flags = {
            component.name: component.validate() for component in self.all_components
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


@dataclass
class GLDSDataSystem(TemplateDataSystem):

    base: BaseDataSystem


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
class TemplateDataset(abc.ABC):
    @property
    @abc.abstractproperty
    def expected_sample_class(self):
        return self.expected_sample_class

    @property
    def samples(self):
        return self.base.samples

    @abc.abstractmethod
    def validate(self):
        ...

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
        sample.dataset = self


#########################################################################
# SAMPLES
#########################################################################
@dataclass
class BaseSample:
    """ Abstract class for samples """

    name: str
    # id: str = field(default_factory=get_id)
    dataset: Union[None, BaseDataset] = field(default=None)


class TemplateSample(abc.ABC):
    @abc.abstractmethod
    def validate(self):
        ...

    # used properties to alias base attributes as needed
    @property
    def name(self):
        return self.base.name

    """ This should be in the mixin
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
        sample.dataset = self
    """


############################################################################################
# BULK RNASEQ SPECIFIC
############################################################################################
class ASSAY(enum.Enum):
    BulkRNASeq = 1


@dataclass
class BulkRNASeqSample(TemplateSample, CanAttachComponents):
    """ Abstract class for samples """

    # composition for all samples
    base: BaseSample
    assay_type: str = ASSAY.BulkRNASeq.name

    # used for paired end
    rawForwardReads: Union[EmptyComponent, ReadsComponent] = field(
        default_factory=EmptyComponent
    )
    rawReverseReads: Union[EmptyComponent, ReadsComponent] = field(
        default_factory=EmptyComponent
    )
    trimForwardReads: Union[EmptyComponent, ReadsComponent] = field(
        default_factory=EmptyComponent
    )
    trimReverseReads: Union[EmptyComponent, ReadsComponent] = field(
        default_factory=EmptyComponent
    )

    # used for single end
    rawReads: Union[EmptyComponent, ReadsComponent] = field(
        default_factory=EmptyComponent
    )
    TrimReads: Union[EmptyComponent, ReadsComponent] = field(
        default_factory=EmptyComponent
    )

    def __post_init__(self):
        pass

    def validate(self):
        strict_type_checks(self)
        # additional checks advised


@dataclass
class BulkRNASeqDataset(TemplateDataset, CanAttachComponents):

    base: BaseDataset
    assay_type: str = ASSAY.BulkRNASeq.name
    expected_sample_class = BulkRNASeqSample

    def validate(self):
        strict_type_checks(self, exceptions=["samples"])


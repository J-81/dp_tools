from ast import Assert
from dataclasses import dataclass, field
import abc
import datetime
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
            self.md5sum = f"DUMMY:{self._compute_md5sum(str(self.path).encode())}"
        else:
            self.md5sum = self._compute_md5sum(self.path.open("rb").read())

        # finally enforce types check
        strict_type_checks(self)

    def _compute_md5sum(self, contents):
        return hashlib.md5(contents).hexdigest()


#########################################################################
# COMPONENTS
#########################################################################
# mixin for attaching samples or datasets or datasystems
@dataclass
class CanAttachEntity:
    entities: dict = field(default_factory=dict)

    def attach_entity(self, entity):
        self.entities[entity.name] = entity


@dataclass
class BaseComponent(CanAttachEntity):
    """ Class for keeping track of abstract components like Reads, Alignments, and Counts """

    # id: str = field(default_factory=get_id)
    description: str = field(default="Null")
    created: datetime.datetime = field(default_factory=datetime.datetime.now)


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
    # id: str = field(default_factory=get_id)
    datasets: Dict[str, "BaseDataset"] = field(default_factory=dict, repr=False)

    # this method should valid the values
    # note: python type checking isn't strictly enforced
    #       So those validations should be performed here
    @abc.abstractmethod
    def __post_init__(self):
        pass

    def attach_dataset(self, dataset: "BaseDataset"):
        if dataset.name in self.datasets.keys():
            log.warning(f"Overwriting pre-existing dataset: {dataset.name}")
        self.datasets[dataset.name] = dataset

    def validate(self):
        # datasystem validations
        ...
        # dataset validations
        [dataset.validate() for dataset in self.datasets.values()]

        # sample validations
        for dataset in self.datasets.values():
            [s.validate() for s in dataset.samples.values()]


@dataclass
class GLDSDataSystem(BaseDataSystem):

    # this method should valid the values
    # note: python type checking isn't strictly enforced
    #       So those validations should be performed here
    def __post_init__(self):
        pass


#########################################################################
# DATASET
#########################################################################
@dataclass
class BaseDataset(abc.ABC):
    """ Abstract class for a dataset
    NOTE: Within this system: 
        a dataset is associated with a specific assay
        a dataSYSTEM is the container for multiple datasets
    """

    # id: str = field(default_factory=get_id)
    name: str = field(default="orphan")
    samples: Dict[str, "BaseSample"] = field(default_factory=dict, repr=False)

    # this method should valid the values
    # note: python type checking isn't strictly enforced
    #       So those validations should be performed here
    @abc.abstractmethod
    def __post_init__(self):
        pass

    @abc.abstractmethod
    def attach_sample(self, sample):
        if sample.name in self.samples.keys():
            log.warning(f"Overwriting pre-existing sample: {sample.name}")
        # attach sample to dataset
        self.samples[sample.name] = sample
        # attach dataset to sample
        sample.dataset = self

    @abc.abstractmethod
    def validate(self):
        strict_type_checks(self, exceptions=["samples"])

        # validate samples
        # for sample in self.samples:
        #    assert isinstance(sample, Sample)


@dataclass
class BulkRNASeqDataset(BaseDataset):

    # sample_class: str = "BulkRNASeqSample"

    def __post_init__(self):
        # append assay type to generate name
        self.name = f"{self.name}:BulkRNASeq"

    def attach_sample(self, sample: "BulkRNASeqSample"):
        try:
            sample.validate()
            super().attach_sample(sample)
        except (AssertionError, AttributeError) as e:
            raise TypeError(
                f"Could not attach {sample} after attempting '{sample}.validate()'"
            )

    def validate(self):
        strict_type_checks(self, exceptions=["samples"])

        # validate samples
        for sample_name, sample in self.samples.items():
            assert sample_name == sample.name
            assert isinstance(sample, BulkRNASeqSample)


#########################################################################
# SAMPLES
#########################################################################
@dataclass
class BaseSample:
    """ Abstract class for samples """

    name: str
    # id: str = field(default_factory=get_id)
    dataset: Union[None, BaseDataset] = field(default=None)

    def validate(self):
        strict_type_checks(self)


# a mixin class
class CanAttachComponents:
    def attach_component(self, component, attr):
        # ensure this is an expected component
        assert hasattr(
            self, attr
        ), "Sample does not specify component with slot: {attr}"
        # warn if component already attached
        if getattr(self, attr):
            log.warn(f"Overwriting existing component in slot: {attr}")

        # attach component
        setattr(self, attr, component)

        # associate component with sample
        component.attach_entity(self)


@dataclass
class BulkRNASeqSample(CanAttachComponents):
    """ Abstract class for samples """

    # composition for all samples
    base: BaseSample

    # used for paired end
    rawForwardReads: Union[None, ReadsComponent] = field(default=None)
    rawReverseReads: Union[None, ReadsComponent] = field(default=None)
    trimForwardReads: Union[None, ReadsComponent] = field(default=None)
    trimReverseReads: Union[None, ReadsComponent] = field(default=None)

    # used for single end
    rawReads: Union[None, ReadsComponent] = field(default=None)
    TrimReads: Union[None, ReadsComponent] = field(default=None)

    def __post_init__(self):
        pass

    def validate(self):
        strict_type_checks(self, except_nones=["dataset"])
        # additional checks advised

    # used properties to alias base attributes as needed
    @property
    def name(self):
        return self.base.name

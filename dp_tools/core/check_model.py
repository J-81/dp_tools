#########################################################################
# FLAG (A validation report message)
#########################################################################
import abc
from dataclasses import dataclass, field
import enum
import os
from pathlib import Path
import re
import traceback
from typing import (
    Callable,
    ClassVar,
    Dict,
    List,
    Optional,
    Set,
    Tuple,
    Union,
)
import logging
from dp_tools.core.entity_model import (
    TemplateComponent,
    TemplateDataset,
    TemplateSample,
)
import pandas as pd
import pkg_resources
from schema import Schema, SchemaMissingKeyError, And, Or
import yaml

log = logging.getLogger(__name__)

from dp_tools.core.model_commons import strict_type_checks

ALLOWED_DEV_EXCEPTIONS = (
    Exception  # Hooking into this with monkeypatch can be handy for testing
)


class FlagCode(enum.Enum):
    """Maps a flag code to a severity level."""

    DEV_HANDLED = 90  # should be used when an error occurs due to developer related mistakes, rather than from the correctly loaded data
    DEV_UNHANDLED = 91  # should never be used directly, instead this should catch unhandled exceptions in the validate function
    HALT1 = 81
    HALT2 = 82
    HALT3 = 83
    HALT4 = 84
    HALT5 = 85
    RED1 = 51
    RED2 = 52
    RED3 = 53
    RED4 = 54
    RED5 = 55
    YELLOW1 = 31
    YELLOW2 = 32
    YELLOW3 = 33
    YELLOW4 = 34
    YELLOW5 = 35
    GREEN = 20
    INFO = 10
    SKIPPED = (
        1  # should never be used directly, instead is the flag code for skipping checks
    )
    DRY_RUN = 0  # should never be used directly, instead is the flag code for validation dry runs

    # allow comparing flag codes
    # used to determine if a code is of higher severity
    def __ge__(self, other):
        return self.value >= other.value

    def __le__(self, other):
        return self.value <= other.value

    def __gt__(self, other):
        return self.value > other.value

    def __lt__(self, other):
        return self.value < other.value


@dataclass
class Check(abc.ABC):

    description: ClassVar[str]
    id: ClassVar[str] = field(init=False)  # set as class name
    flag_desc: ClassVar[Dict[int, str]]
    # flags associated with this check
    flags: List["Flag"] = field(default_factory=list)
    config: ClassVar[Dict] = dict()  # set as empty dict if a check doesn't need config
    allChecks: ClassVar[
        List["Check"]
    ] = list()  # note this is a class attribute list, tracking all checks
    _dry_run: ClassVar[bool] = field(
        default=False, repr=False
    )  # controls whether a validation is performed as per dry run
    skip_this: ClassVar[bool] = field(
        default=False, repr=False
    )  # controls whether a validation is performed as per explicit skipping

    def __post_init__(self, skip: bool = None):
        # format description with config
        # All templated fields MUST be in config
        # All items in config are NOT required
        # (however, they should be supplied to make the description complete)
        self._proto_description = self.description
        # assert that this is a string and not a tuple on accident
        assert not isinstance(
            self._proto_description, tuple
        ), "Check description must be a string, check that you did not add commas to a multi-line str tuple"
        assert isinstance(self._proto_description, str)
        self.description = self.description.format(**self.config)

        # set check id as class name
        self.id = self.__class__.__name__

        # ensure checkID is valid
        assert re.match(
            r"^(COMPONENT|SAMPLE|DATASET)_\D*_\d\d\d\d$", self.id
        ), "Valid Check ID format: SCOPE_WORD_dddd where d is any digit. SCOPE must be either 'COMPONENT','SAMPLE' or 'DATASET'"

        # ensure check is unique
        assert self.id not in [
            check.id for check in self.allChecks
        ], "CheckID already used, try another ID"

        # add to tracked allChecks
        self.allChecks.append(self)

        # set skip if added
        if skip != None:
            self.__class__.skip_this = skip

        # add flag description for unhandled developer exception flags
        self.flag_desc[
            FlagCode.DEV_UNHANDLED
        ] = "An unhandled exception occured in {function} with args '{args}' and kwargs: '{kwargs}': {e}.\n{full_traceback}"

        # add flag description for handling DRYRUN flags
        self.flag_desc[FlagCode.DRY_RUN] = "DRY RUN REQUESTED: check not run"

        # add flag description for handling SKIPPED flags
        self.flag_desc[FlagCode.SKIPPED] = "CHECK SKIP REQUESTED: check not run"

    @abc.abstractmethod
    def validate_func(self) -> "Flag":
        raise NotImplementedError

    def validate(self, *args, request_skip: bool = False, **kwargs):
        # add skipped route that bails out
        # NOTE: this intentionally takes place before the dry run condition is checked
        #       this allows checking if your config will perform as expected, including skips
        if request_skip:
            return Flag(
                codes=FlagCode.SKIPPED,
                message_args={},
                check=self,
            )

        # add dry run route that bails out
        if self.__class__._dry_run:
            return Flag(
                codes=FlagCode.DRY_RUN,
                message_args={},
                check=self,
            )

        try:
            return self.validate_func(*args, **kwargs)
        except ALLOWED_DEV_EXCEPTIONS as e:
            log.critical(
                f"Developer exception raised during a validation function for check: {self}"
            )
            return Flag(
                codes=FlagCode.DEV_UNHANDLED,
                message_args={
                    "function": self.validate_func,
                    "e": e,
                    "full_traceback": traceback.format_exc(),
                    "args": args,
                    "kwargs": kwargs,
                },
                check=self,
            )

    def copy_with_new_config(self, id: str, config: dict) -> "Check":
        return Check(
            id=id,
            config=config,
            description=self._proto_description,
            validate_func=self.validate_func,
            flag_desc=self.flag_desc,
        )

    def __str__(self):
        return f"""ID: {self.id} # of Flags Generated: {len(self.flags)}
        description: {self.description}
        """


@dataclass
class Flag:

    codes: set[FlagCode]
    check: Check
    message_args: dict = field(default_factory=dict)
    message: str = field(init=False)

    # must ensure all flags are correct before continuing
    def __post_init__(self):

        # code related init
        # covert single codes to list if not already a list
        if isinstance(self.codes, FlagCode):
            self.codes = {self.codes}

        # remove FlagCode.GREEN if included with other flags
        # that flag should never be reported if higher severity flags are reported
        if FlagCode.GREEN in self.codes and max(self.codes) != FlagCode.GREEN:
            self.codes.remove(FlagCode.GREEN)

        # convert codes into sorted list in descending severity order
        self.codes = list(self.codes)
        self.codes.sort(reverse=True)

        # set message
        self.message = None
        for code in self.codes:
            if not self.message:
                self.message = ""
            else:
                self.message += ":::"
            # retrieve proto message
            pre_format_msg = self.check.flag_desc[code]

            # populate message with message_args
            self.message += pre_format_msg.format(**self.message_args)

        # check types before final addition and init
        strict_type_checks(
            self, exceptions=["codes"]
        )  # codes is a list of objects, not currently covered by this check

        # add to check flags records
        self.check.flags.append(self)

    @property
    def maxCode(self):
        """Returns most severe FlagCode rather than all codes for backwards compatibility"""
        return max(self.codes)

    def __str__(self):
        return f"Flag (checkID: {self.check.id}, FlagCodes: {self.codes}, message: {self.message})"


class VVProtocol(abc.ABC):
    """A abstract protocol for VV"""

    def __init__(
        self,
        dataset,
        config: Union[Tuple[str, str], Path],
        dry_run: bool = False,
        skip_these_checks: set = None,
        stage_names=None,  # TODO: fix typehint here
        protocol_name: str = "DEFAULT",
    ):

        # logging
        log.info(
            f"Initiating V&V protocol with config: '{config}' and '{protocol_name}'"
        )

        # parse config for values
        # TODO: raise exceptions/warnings when both CLI and config args are double supplied
        if validation_config := self.load_config(config, protocol_name):
            try:
                # retrieve from config
                # at this point, this is a string
                config_stage_names: Union[str, Set[str]] = validation_config["stage"]
            except KeyError:
                raise KeyError(
                    f"Config file did not specify a valid stage: valid={[i for i in self.STAGES]}, supplied:{validation_config['stage']}"
                )

            if skip_these_checks == None:
                skip_these_checks = validation_config.get("skip these checks", None)

        # set as final stage_names if not overridden
        # TODO: invert this for clarity
        if not stage_names:
            stage_names = config_stage_names

        # stage can be either a list of stages or a single stage
        # for a single stage, all stages preceeding should be included
        if isinstance(stage_names, str):
            self._stage_arg = self.STAGES[stage_names]
            self._stages_loaded = {self._stage_arg}.union(
                self.STAGES.get_all_preceeding(self._stage_arg)
            )
        elif isinstance(stage_names, set):
            self._stage_arg = {self.STAGES[stage_name] for stage_name in stage_names}
            # validate contents
            for sub_stage in self._stage_arg:
                if not isinstance(sub_stage, self.STAGES):
                    raise ValueError(f"All requested stages must be a {self.STAGES}")
            # set as loaded stages
            self._stages_loaded = self._stage_arg
        else:
            raise ValueError(
                "'stage' must be supplied either directly or as part of configuration"
            )

        log.info(f"Loaded validation stages: {self._stages_loaded}")

        # assign to empty set if no skips are specified
        if skip_these_checks == None:
            self.skip_these_checks = set()
        else:
            self.skip_these_checks = skip_these_checks
        # check on init that dataset is of the right type
        assert isinstance(
            dataset, self.expected_dataset_class
        ), f"dataset MUST be of type {self.expected_dataset_class} for this protocol"
        self.dataset = dataset
        self._flags = {"dataset": dict(), "sample": dict(), "component": dict()}
        if dry_run:
            log.critical("DRY RUNNING: Validation protocol, no checks performed")
            Check._dry_run = True

    def run_check(self, check: Check, *args, **kwargs) -> Flag:
        """Runs check and performs validation if skip is not requested"""
        return check.validate(
            request_skip=(check.id in self.skip_these_checks), *args, **kwargs
        )

    @property
    def flags(self):
        return self._flags

    @abc.abstractproperty
    @property
    def expected_dataset_class(self):
        return self.expected_dataset_class

    @abc.abstractproperty
    @property
    def STAGES(self):
        return self._STAGES

    def validate_all(self):
        self._flags["dataset"].update(self.validate_dataset())
        self._flags["sample"].update(self.validate_samples())
        self._flags["component"].update(self.validate_components())

    @abc.abstractmethod
    def validate_dataset(self) -> Dict[TemplateDataset, List[Flag]]:
        ...

    @abc.abstractmethod
    def validate_samples(self) -> Dict[TemplateSample, List[Flag]]:
        ...

    @abc.abstractmethod
    def validate_components(self) -> Dict[TemplateComponent, List[Flag]]:
        ...

    def load_config(
        self, config: Union[Tuple[str, str], Path], protocol_name: str = "DEFAULT"
    ) -> dict:
        if isinstance(config, tuple):
            conf_validation = yaml.safe_load(
                pkg_resources.resource_string(
                    __name__,
                    os.path.join("..", "config", f"{config[0]}_v{config[1]}.yaml"),
                )
            )["Validation Protocols"][protocol_name]
        elif isinstance(config, Path):
            conf_validation = yaml.safe_load(config.open())["Validation Protocols"][
                protocol_name
            ]

        log.debug("Loaded the following validation config: {conf_validation}")

        # validate with schema
        config_schema = Schema(
            {"dry run": bool, "skip checks": list, "stage": Or(list[str], str)}
        )

        config_schema.validate(conf_validation)

        # convert list to sets
        conf_validation["skip checks"] = set(conf_validation["skip checks"])
        # convert only if this is a list
        conf_validation["stage"] = (
            set(conf_validation["stage"])
            if isinstance(conf_validation["stage"], list)
            else conf_validation["stage"]
        )

        self.config = conf_validation
        return conf_validation

    ###################################################################################################################
    # EXPORT RELATED METHODS (POST VALIDATION)
    ###################################################################################################################
    def flags_to_df(self, schema="default") -> pd.DataFrame:
        log.debug(f"Extracting info into schema {schema} records")
        match schema:
            # default dataframe format
            # entity_index, flag_code, checkID, message
            # tuple[str, str, str, str], int, str, str
            case "default":
                records = list()
                for entity_type in ["dataset", "sample", "component"]:
                    for entity, flags in self.flags[entity_type].items():
                        record_partial = {
                            "entity_index": self._get_entity_index(entity)
                        }
                        for flag in flags:
                            record_full = record_partial | {
                                "flag_code": flag.maxCode.value,
                                "check_id": flag.check.id,
                                "message": flag.message,
                            }
                            records.append(record_full)
            # verbose dataframe format
            # entity_index, flag_enum, flag_code, message, checkID, checkDescription
            # tuple[str, str, str, str], str, int, str, str
            case "verbose":
                records = list()
                for entity_type in ["dataset", "sample", "component"]:
                    for entity, flags in self.flags[entity_type].items():
                        record_partial = {
                            "entity_index": self._get_entity_index(entity)
                        }
                        for flag in flags:
                            record_full = record_partial | {
                                "flag_name": flag.maxCode.name,
                                "flag_code": flag.maxCode.value,
                                "all_flag_codes": flag.codes,
                                "message": flag.message,
                                "check_id": flag.check.id,
                                "check_description": flag.check.description,
                            }
                            records.append(record_full)
            case _:
                raise ValueError(f"Schema: {schema} not defined")
        log.debug("Exporting records to dataframe")
        # set up datafrom with simple index
        df = pd.DataFrame.from_records(data=records).set_index("entity_index")
        # convert into multiIndex
        df.index = pd.MultiIndex.from_tuples(
            df.index, names=["DataSystem", "Dataset", "Sample", "Component"]
        )
        return df

    def _get_entity_index(
        self, entity: Union[TemplateDataset, TemplateSample, TemplateComponent]
    ) -> Tuple[str, str, str, str]:
        match entity:
            case TemplateDataset():
                entity_index = (entity.dataSystem.name, entity.name, "", "")
            case TemplateSample():
                entity_index = (
                    entity.dataset.dataSystem.name,
                    entity.dataset.name,
                    entity.name,
                    "",
                )
            case TemplateComponent(entityType="sample"):
                sample_entities = ",".join([sample for sample in entity.entities])
                attr_index = ",".join([d["attr"] for d in entity.entities.values()])
                entity_index = (
                    entity.dataset.dataSystem.name,
                    entity.dataset.name,
                    sample_entities,
                    attr_index,
                )
            case TemplateComponent(entityType="dataset"):
                attr_index = ",".join([d["attr"] for d in entity.entities.values()])
                entity_index = (
                    entity.dataset.dataSystem.name,
                    entity.dataset.name,
                    "",
                    attr_index,
                )
            case _:
                raise ValueError(f"Unable to extract entity index from : {entity}")
        return entity_index

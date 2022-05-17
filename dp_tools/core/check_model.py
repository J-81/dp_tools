#########################################################################
# FLAG (A validation report message)
#########################################################################
import abc
from collections import defaultdict
from contextlib import contextmanager
from dataclasses import dataclass, field
import enum
import json
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
    TypedDict,
    Union,
)
import logging
#from typing_extensions import NotRequired
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


'''
@dataclass
class Check(abc.ABC):

    USE_SUBCHECK_MONADS: ClassVar[str] = field(
        default=False
    )  # defines whether to use subcheck monads and associated formatting approaches
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
        try:
            self.__class__.description = self.description.format(**self.config)
        except Exception as e:
            log.error(
                f"Exception occured during description formatting: description_template: '{self.description}', configuration: '{self.config}'"
            )
            raise

        # set check id as class name
        self.__class__.id = self.__class__.__name__

        # ensure checkID is valid
        assert re.match(
            r"^(COMPONENT|SAMPLE|DATASET)_\D*_\d\d\d\d$", self.id
        ), "Valid Check ID format: SCOPE_WORD_dddd where d is any digit. SCOPE must be either 'COMPONENT','SAMPLE' or 'DATASET'"

        # ensure check is unique
        assert self.id not in [
            check.id for check in self.allChecks
        ], f"CheckID {self.id} already used, try another ID"

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
        log.info(f"Running check: '{self.id}' with args: {args}, kwargs: {kwargs}")
        # add skipped route that bails out
        # NOTE: this intentionally takes place before the dry run condition is checked
        #       this allows checking if your config will perform as expected, including skips
        if request_skip:
            self.USE_SUBCHECK_MONADS = False  # this reverts back to Flag with codes as FlagCodes rather than a flag_data stream
            log.info(f"  Skipping as per request!")
            return Flag(
                codes=FlagCode.SKIPPED,
                message_args={},
                check=self,
            )

        # add dry run route that bails out
        if self.__class__._dry_run:
            self.USE_SUBCHECK_MONADS = False  # this reverts back to Flag with codes as FlagCodes rather than a flag_data stream
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
            self.USE_SUBCHECK_MONADS = False  # this reverts back to Flag with codes as FlagCodes rather than a flag_data stream
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


class FlagEntry(TypedDict):
    # programmed into validation function or mapping function
    code: FlagCode
    message: str
    # attached by monad
    args: list
    kwargs: dict
    check_function_name: Callable


class Flaggable(object):
    def __init__(self, *args, flag_data: list[FlagEntry] = None, **kwargs):
        self.args = args
        self.kwargs = kwargs
        self.flag_data = flag_data if flag_data is not None else list()
        log.debug(f"Created monad {self}")

    def __repr__(self):
        return f"Flaggable({self.args}, {self.kwargs}, {self.flag_data})"

    def flag_data_to_df(self):
        return pd.DataFrame(self.flag_data)

    def bind(
        self,
        f: Callable,
        skip: bool = False,
        *args,
        **kwargs,
    ) -> "Flaggable":
        log.debug(f"Binding function {f} to monad {self}")
        # remove bind specific kwargs
        combined_args = self.args + args
        combined_kwargs = self.kwargs | kwargs
        assert isinstance(skip, bool), "bind argument 'skip' MUST evaluate to a boolean"
        if skip == True:
            result = {
                "code": FlagCode.SKIPPED,
                "message": f"Skipped as the bind argument 'skip' evaluated to 'True'",
                "args": combined_args,
                "kwargs": combined_kwargs,
                "check_function_name": f.__name__,
            }
            updated_flag_data = self.flag_data + [result]
            return_monad = Flaggable(
                *self.args, **self.kwargs, flag_data=updated_flag_data
            )
            return return_monad
        try:
            # filter to only keywords args that the function uses
            # f_kwargs = list(inspect.signature(f).parameters.keys())
            # f_payload = {k: v for k, v in combined_kwargs.items() if k in f_kwargs}
            # result = f(*combined_args, **f_payload)
            result = f(*combined_args, **combined_kwargs)
            result = result | {
                "args": combined_args,
                "kwargs": combined_kwargs,
                "check_function_name": f.__name__,
            }
            updated_flag_data = self.flag_data + [result]
            return_monad = Flaggable(
                *self.args, **self.kwargs, flag_data=updated_flag_data
            )
            log.debug(f"Returning computed function {f} resulting in {return_monad}")
            return return_monad
        except AssertionError as e:
            raise ValueError(f"Bad return: {result}") from e
        except SystemExit as e:
            failure_status = {
                "trace": traceback.format_exc(),
                "exc": e,
                "args": args,
                "kwargs": kwargs,
            }
            result = {
                "code": FlagCode.DEV_UNHANDLED,
                "message": f"Raised an unhandled exception: {failure_status}",
                "args": combined_args,
                "kwargs": combined_kwargs,
                "check_function_name": f.__name__,
            }
            updated_flag_data = self.flag_data + [result]
            return_monad = Flaggable(
                *self.args, **self.kwargs, flag_data=updated_flag_data
            )
            log.debug(f"Returning computed function {f} resulting in {return_monad}")
            return return_monad

    def export_flag_data_to_Flag(self, check: Check) -> "Flag":
        """Exports flag data to generate return a Flag

        :return: Flag object compatible with the check_model
        :rtype: Flag
        """

        return Flag(
            codes=self.flag_data, check=check, message_args={}  # unused in Flag
        )


@dataclass
class Flag:

    codes: Union[
        set[FlagCode], list[FlagEntry]
    ]  # allow use of flag data from Flaggable monad
    check: Check
    message_args: dict = field(default_factory=dict)
    message: str = field(init=False, default=None)

    # must ensure all flags are correct before continuing
    def __post_init__(self):
        if self.check.USE_SUBCHECK_MONADS:
            # extract message
            self.message = json.dumps(
                {
                    e["check_function_name"]: {
                        "code": e["code"].name,
                        "message": e["message"],
                    }
                    for e in self.codes
                }
            )

            # extract flag codes
            self.codes = {e["code"] for e in self.codes.copy()}

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

        # set message, if not created by flag monad
        if self.message is None:
            self.message = "{"
            for code in self.codes:
                # retrieve proto message
                pre_format_msg = self.check.flag_desc[code]
                # populate message with message_args
                add_string = f"{pre_format_msg.format(**self.message_args)}"
                self.message += f"'{code.name}':'{add_string}',"
            # the -1 removes the last comma
            self.message = self.message[:-1] + "}"  # to close the json like message

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
        self._flags = {
            "dataset": defaultdict(list),
            "sample": defaultdict(list),
            "component": defaultdict(list),
        }
        if dry_run:
            log.critical("DRY RUNNING: Validation protocol, no checks performed")
            Check._dry_run = True

        self.check_cache = dict()

    def validate(self, level: str, **entity) -> "VVProtocol":
        """Bind the entities to validate.
        *args and **kwargs are passed to all included validation checks

        :return: The protocol after binding entities (for chain methods)
        :rtype: VVProtocol
        """
        # TODO: add logging re: overwriting bound entities
        assert len(entity) == 1, "Only one bound entity is allowed"
        self.entity = entity
        self.level = level
        return self

    def include(self, check: Check) -> "VVProtocol":
        """Bind a check to internally run and append to flags

        :param check: A check that returns a flag when run
        :type check: Check
        :return: The protocol after running the check and appending the returned flag to flags
        :rtype: VVProtocol
        """

        # initiate check
        check_to_run: Check = check()
        [entity] = (entity for entity in self.entity.values())
        self.flags[self.level][entity] = (
            self.flags[self.level].get(entity)
            if self.flags[self.level].get(entity)
            else list()
        )
        self.flags[self.level][entity].append(
            check_to_run.validate(
                request_skip=(check_to_run.id in self.skip_these_checks), **self.entity
            )
        )
        return self

    def run_check(self, check: type[Check], *args, **kwargs) -> Flag:
        """Runs check and performs validation if skip is not requested"""
        return check().validate(
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
        self.validation_procedure()
        self._flags["component"].update(self.validate_components())

    @abc.abstractmethod
    def validation_procedure(self):
        ...

    def safe_init_check(self, check: type[Check]):
        if cached_check := self.check_cache.get(check.__name__):
            return cached_check
        else:
            self.check_cache[check.__name__] = check()
            return self.check_cache[check.__name__]

    def batch_check(
        self,
        level: str,
        payloads: list[dict],
        run: list[type[Check]],
        payload_as_key: bool = True,
    ):
        log.info(f"Starting batch check at '{level}' level.")
        for check in run:
            init_check = self.safe_init_check(check)
            for payload in payloads:
                log.info(f"Running {check.__name__} with payload: {payload}")
                flag = init_check.validate(**payload)
                if payload_as_key:
                    [payload_value] = (v for v in payload.values())
                    self._flags[level][payload_value].append(flag)

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
                                "flag_name": flag.maxCode.name,
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
'''

########################################################################
########################################################################
########################################################################
########################################################################
class FlagEntry(TypedDict):
    code: FlagCode
    message: str

class FlagEntryWithOutliers(FlagEntry):
    # TODO: make typehint
    outliers: dict[str, dict[str, dict[str, str]]]


class ValidationProtocol:
    class Component:
        def __init__(
            self,
            name: str,
            parent: "ValidationProtocol.Component",
            skip: bool = False,
            description: str = "",
        ):
            log.info(f"Creating new entity in validation protocol: {name}")
            self.flags: list = list()
            self.parent = parent
            self.description = description
            self.name = name
            # only root entity should have no parent
            if self.parent:
                self.parent.children.append(self)
            self.children: list["ValidationProtocol.Component"] = list()
            self.skip = skip

        def __repr__(self):
            if self.name == "ROOT":
                return f"Component(name={self.name}, parent=!!NONE THIS IS THE ROOT COMPONENT!!)"
            return f"Component(name={self.name}, parent={self.parent.name})"

        @property
        def ancestor_line(self):
            if self.parent is not None:
                return self.parent.ancestor_line + [self.name]
            else:
                return [self.name]

        def ancestry_is_in(self, other_list):
            """Returns True if the component or any ancestor of the component is
            in the 'other_list'
            """
            return any([(name in other_list) for name in self.ancestor_line])

    class QueuedCheck(TypedDict):
        check_fcn: Callable
        description: str
        payload: dict
        config: dict
        component: "ValidationProtocol.Component"
        to_run: bool

    class _ALL_COMPONENTS(list):
        """Dummy list that works as a default whitelist of all components"""

        def __init__(self):
            ...

        def __contains__(self, _):
            return True

    def __init__(self, skip_components: list = None, run_components: list = None):
        if run_components is not None:
            self.run_components = run_components
        else:
            self.run_components = self._ALL_COMPONENTS()

        if skip_components is None:
            skip_components = list()

        self.skip_components = skip_components
        # init
        # self._stage = defaultdict(list)
        self._payloads: list = list()
        # TODO: typehint dicts here
        self.outliers: dict[str, dict[str, dict[str, dict[str, str]]]] = defaultdict(
            lambda: defaultdict(lambda: defaultdict(dict))
        )

        def nested_defaultdict():
            return defaultdict(nested_defaultdict)

        self.results: dict = defaultdict(nested_defaultdict)
        # typehint exception: only the root component breaks the MUST have a parent rule
        self._root_component = ValidationProtocol.Component(name="ROOT", parent=None)  # type: ignore
        self.cur_component = self._root_component

        # will hold all pending checks
        self._check_queue: list[ValidationProtocol.QueuedCheck] = list()

    @contextmanager
    def component_start(self, name: str, description: str, skip: bool = False):
        try:
            # set on protocol init
            # set current stage as parent
            new_component = ValidationProtocol.Component(
                name=name, parent=self.cur_component, description=description
            )
            self.cur_component = new_component

            # determine if skipping needs to be set for the new component
            if any(
                [
                    new_component.ancestry_is_in(self.skip_components),
                    not new_component.ancestry_is_in(self.run_components),
                    skip,
                ]
            ):
                new_component.skip = True
            yield
        finally:
            self.cur_component = self.cur_component.parent

    @staticmethod
    def _eval_payload_callables(payload: dict) -> dict:
        evaluated_payload = dict()
        for k, v in payload.items():
            if callable(v):
                evaluated_payload[k] = v()
            else:
                evaluated_payload[k] = v
        return evaluated_payload

    @contextmanager
    def payload(self, payloads: list[dict]):
        try:
            self._payloads = payloads
            yield self.add
        finally:
            # clear payloads on exit
            self._payloads = list()

    def add(
        self,
        fcn: Callable[..., FlagEntry],
        payloads: dict = None,
        skip: bool = False,
        config: dict = None,
        description: str = None,
    ):
        # override payload with one supplied directly to run
        if payloads:
            self._payloads = [payloads]

        # coerce config variable if not provided to an empty dict
        if config is None:
            config = dict()

        for payload in self._payloads:
            self._check_queue.append(
                {
                    "check_fcn": fcn,
                    "description": description
                    if description is not None
                    else fcn.__name__,
                    "payload": payload,
                    "config": config,
                    "component": self.cur_component,
                    # don't run if either this add call specifies skip
                    # or if the component is being skipped
                    "to_run": not any([skip, self.cur_component.skip]),
                }
            )
        return self

    def run(self):
        """Runs all queue checks"""
        for queued in self._check_queue:
            fcn = queued["check_fcn"]
            fcn_name = fcn.__name__

            if not queued["to_run"]:
                # case: skipping compute
                packed_result = {
                    "code": FlagCode.SKIPPED,
                    "message": "Skipped by protocol",
                    "function": fcn_name,
                    "description": queued["description"],
                    # NOTE: This is intentionally the unevaluated payload
                    #     This assumes skipped checks may contain inputs that fail on evaluation
                    #     One example being data models that aren't fully loaded for skipped checks
                    "kwargs": queued["payload"] | queued["config"],
                    "config": queued["config"],
                }
            else:
                # case: to_compute is true
                # try computing, works
                payload = self._eval_payload_callables(queued["payload"])
                payload_and_config = payload | queued["config"]
                try:
                    result = fcn(**payload_and_config)

                    # peel off outlier data and keep track
                    # using current component name as the top level key
                    if fcn_outliers := result.pop("outliers", None):
                        self.outliers[queued["component"].name] = (
                            self.outliers[queued["component"].name] | fcn_outliers
                        )

                    # note: we want to make sure to override description if the result gives one;
                    # however it is optional and the fallback is the function name
                    packed_result = (
                        {"description": queued["description"]}
                        | result
                        | {
                            "kwargs": payload_and_config,
                            "config": queued["config"],
                            "function": fcn_name,
                        }
                    )
                # except: computing failed
                # on exception, capture exception
                except SystemExit as e:
                    packed_result = {
                        "code": FlagCode.DEV_UNHANDLED,
                        "message": e,
                        "function": fcn_name,
                        "description": queued["description"],
                        "kwargs": payload_and_config,
                        "config": queued["config"],
                    }
            # add result (including skip flag) to component
            queued["component"].flags.append(packed_result)

    ##############################################
    ### METHODS FOR EXPORTING FLAG DATA
    ##############################################

    @staticmethod
    def get_report_data(component):
        yield {tuple(component.ancestor_line): component.flags}
        for child in component.children:
            yield from ValidationProtocol.get_report_data(child)

    def report(
        self,
        nested: bool = True,
        include_flags: bool = True,
        include_skipped: bool = True,
    ) -> dict:
        """Returns a report of all components and their associated flags

        :param nested: If True, the resultant dictionary is nested analogous to component tree
        :param include_flags: If True, the resultant dictionary includes flags
        :param include_skipped: If True, the resultant dictionary includes skipped flags
        :return: Dictionary of all components by name and their flags
        :rtype: dict
        """
        self.results = dict(self.results)
        data = list(self.get_report_data(self._root_component))
        df_data: list[dict] = list()
        for component in data:
            for index, flags in component.items():
                for flag in flags:
                    entry = {
                        "index": index,
                    } | flag
                    df_data.append(entry)

        df = pd.DataFrame(df_data).set_index("index")

        # filtering as requested in args
        if not include_skipped:
            df = df.loc[df["code"] != FlagCode.SKIPPED]

        def default_to_regular(d):
            if isinstance(d, defaultdict):
                d = {k: default_to_regular(v) for k, v in d.items()}
            return d

        return {
            "flag_table": df,
            "outliers": pd.DataFrame.from_dict(
                default_to_regular(self.outliers), orient="index"
            ).T,
        }


# TODO: Part of an alternative framework not fully implemented
"""
@dataclass
class SingleFlag:
    # expected as check function return
    code: FlagCode
    message: str
    # attached by protocol
    check_function: Callable
    description: str

class Entity:
    def __init__(self, data: object, parent: "Entity"):
        log.info(f"Creating new entity in validation protocol: {data}")
        self.data = data
        self.flags: List[SingleFlag] = list()
        self.parent = parent
        # only root entity should have no parent
        if self.parent:
            self.parent.children.append(self)
        self.children: list["Entity"] = list()


class Protocol(abc.ABC):
    def __init__(self):
        self.cursor = None
        self.skipping = False

    # Accepts single entity as a kwarg to ensure pass down to checks as kwargs
    def start_component(self, **data):
        log.info(f"Starting validation, setting root entity as {data}")
        entity = Entity(data=data, parent=self.cursor)
        # set root entity if this is the first component
        if self.cursor is None:
            self.root = entity
        self.cursor = entity

        return self

    def run(
        self,
        func_check: Callable,
        iff_parent_passing: bool = False,
    ):
        if iff_parent_passing and not all(self.cursor.parent.flags):
            log.warn(f"Parent flags were not all passing, skipping the next block")
            self.cursor.flags.append("SKIPPED!")
            return self

        log.info(f"Running {func_check} with data: {self.cursor.data}")
        code, message = func_check(**self.cursor.data)

        flag = SingleFlag(code=code, message=message, check_function=func_check)
        self.cursor.flags.append(flag)

        return self

    def end_component(self):
        log.info(
            f"Ending {self.cursor} validation, moving back to {self.cursor.parent}"
        )
        self.cursor = self.cursor.parent

        return self
"""

#########################################################################
# FLAG (A validation report message)
#########################################################################
import abc
from collections import Counter, defaultdict
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
    """Maps a flag code to a severity level.

    Note:
        Severity of flag codes are as follows:

        DEV_UNHANDLED > HALT > RED > YELLOW > GREEN > INFO > SKIPPED

    Example:
        Comparing flags will compare severity levels::

            FlagCode.RED > FlagCode.GREEN
            # should return True

    """

    DEV_UNHANDLED = 90
    """Never used directly, instead used when an unhandled exceptions is raised in the check function"""
    HALT = 80
    """Denotes an issue problematic enough to merit halting an ongoing process immediately. Note: Actually issuing any halt commands is outside the scope of validation protocol itself."""
    RED = 50
    """A severe issue. Likely requiring further manual investigation."""
    YELLOW = 30
    """A moderate/minor issue. May require further manual investigation especially if concentrated around a specific component."""
    GREEN = 20
    """Indicates no issues found during the check."""
    INFO = 10
    """Denotes the results of a check without an issue catching mechanism. May be used to report information in the flag table."""
    SKIPPED = 1
    """Never used directly, instead is the flag code for skipped checks"""

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


########################################################################
########################################################################
########################################################################
########################################################################
class FlagEntry(TypedDict):
    """A dictionary reporting results of a check."""

    code: FlagCode
    """ Assessment of the check's success/failure """
    message: str
    """ A check specific message describing why the code was issued """


class FlagEntryWithOutliers(FlagEntry):
    """Extension of FlagEntry that includes outliers information."""

    # TODO: make typehint
    outliers: dict[str, dict[str, dict[str, str]]]
    """ A dictionary describing outliers found during the check """


class ValidationProtocol:
    """A validation protocol.

    Validation protocols are used to connect check function, function configuration
    and inputs to form a check queue.

    The check queue can then be printed and/or run.

    Finally, the protocol object includes reporting functions
    to present the results of those checks.

    Example:
        An entire sample validation protocol run::

            from dp_tools.core.check_model import ValidationProtocol, FlagCode

            car = {
                "wheels": [
                    {"isFlat": False},
                    {"isFlat": False},
                    {"isFlat": True},  # uh oh
                    {"isFlat": False},
                ],
                "engineOkay": True,
                "gasTank": 100,
                "isElectric": False,
                # doesn't apply so this entry isn't included
                # "charge":0
            }


            def check_if_wheel_flat(wheel: dict):
                if not wheel["isFlat"]:
                    code = FlagCode.GREEN
                    message = "Engine looks good"
                else:
                    code = FlagCode.HALT1
                    message = "Engine needs work!"
                return {"code": code, "message": message}


            def check_engine_okay(car: dict):
                if car["engineOkay"]:
                    code = FlagCode.GREEN
                    message = "Engine looks good"
                else:
                    code = FlagCode.HALT1
                    message = "Engine needs work!"
                return {"code": code, "message": message}


            def check_gas(cur_gas: float, minimum_gas: float):
                if cur_gas >= minimum_gas:
                    code = FlagCode.GREEN
                    message = f"Gas tank is at {cur_gas}. Which is above {minimum_gas}"
                else:
                    code = FlagCode.HALT1
                    message = (
                        f"Gas tank needs a fill up! Current: {cur_gas}. Minimum: {minimum_gas}"
                    )
                return {"code": code, "message": message}


            def check_charge(cur_charge: float, minimum_charge: float):
                if cur_charge >= minimum_charge:
                    code = FlagCode.GREEN
                    message = "Charge looks good"
                else:
                    code = FlagCode.HALT1
                    message = "Needs a charge!"
                return {"code": code, "message": message}


            vp = ValidationProtocol()
            with vp.component_start(
                name="PreDrive Check",
                description="Make sure the car is ready for the trip",
            ):

                with vp.component_start(
                    name="Engine",
                    description="Make sure the engine is running fine",
                ):
                    with vp.payload(payloads=[{"car": car}]):
                        vp.add(check_engine_okay)

                with vp.component_start(
                    name="Tires",
                    description="Make sure the tires are road ready",
                ):
                    with vp.payload(
                        payloads=[
                            {"wheel": car["wheels"][0]},
                            {"wheel": car["wheels"][1]},
                            {"wheel": car["wheels"][2]},
                            {"wheel": car["wheels"][3]},
                        ]
                    ):
                        vp.add(check_if_wheel_flat)

                with vp.component_start(
                    name="Fuel",
                    description="Check gas or charge is sufficent",
                ):
                    # NOTE: lambda is used in payload to defer evaluation conditioned
                    #   on whether the check is run or skipped.
                    with vp.payload(payloads=[{"cur_gas": lambda: car["gasTank"]}]):
                        vp.add(check_gas, config={"minimum_gas": 95}, skip=(car["isElectric"]))

                    with vp.payload(payloads=[{"cur_charge": lambda: car["charge"]}]):
                        vp.add(
                            check_charge,
                            config={"minimum_charge": 95},
                            skip=(not car["isElectric"]),
                        )

                # Now printing the queued checks
                print(vp.queued_checks())

                print("********************************")
                # Running the checks
                vp.run()

                # And showing the results
                print(vp.report()["flag_table"])

                print("Individual Flag Row")
                print(vp.report()["flag_table"].iloc[5].to_dict())

    Args:
        skip_components (list, optional): Component names to skip. Defaults to None.
        run_components (list, optional): Components name to run. Defaults to None.
            If no components are specified, all components will run (except those explicitly
            included in the 'skip_components')
    """

    class _Component:
        def __init__(
            self,
            name: str,
            parent: "ValidationProtocol._Component",
            skip: bool = False,
            skip_children: bool = False,
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
            self.children: list["ValidationProtocol._Component"] = list()

            self.skip = skip
            """ Denotes if this components checks will be skipped """
            self.skip_children = skip_children
            """ Denotes if the components children should be skipped, propogates to their children as well """

        def __repr__(self):
            if self.name == "ROOT":
                return f"Component(name={self.name}, parent=!!NONE THIS IS THE ROOT COMPONENT!!)"
            return f"Component(name={self.name}, parent={self.parent.name}, skip={self.skip}, skip_children={self.skip_children})"

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

    class _QueuedCheck(TypedDict):
        """A queued check including checks that will be skipped"""

        check_fcn: Union[str, Callable[..., FlagEntry]]
        """ A callable function that returns a flag entry or a string placeholder"""

        fcn_name: str
        """ Denotes the function or placeholder name """

        description: str
        """ A user friendly description of the specific check, useful when a general check function where the function name alone is too non-specific"""

        payload: dict
        """ The keyword arguments that will be supplied to the function if it is run. Considered 'dynamic' compared to the config (e.g. the specifc log file path for a validation protocol execution is best defined as a payload) """

        config: dict
        """ Additional keyword arguments that will be supplied to the function. Considered 'static' compared to payload (e.g. 'number of lines to check' is best defined as check configuration)"""

        component: "ValidationProtocol._Component"
        """ The component that results of the check will be attached to for reporting purposes """

        to_run: bool
        """ Defines whether the check will be executed on protocol run (and whether the payload will be evaluated in the case of lambda style payloads). """

    class _ALL_COMPONENTS(list):
        """Dummy list that works as a default whitelist of all components"""

        def __init__(self):
            ...

        def __contains__(self, _):
            return True

        def __repr__(self):
            return "ALL_COMPONENTS"

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
        self._root_component = ValidationProtocol._Component(name="ROOT", parent=None)  # type: ignore
        self.cur_component = self._root_component

        # will hold all pending checks
        self._check_queue: list[ValidationProtocol._QueuedCheck] = list()

    ##############################################
    ### METHODS FOR PLANNING VALIDATION CHECKS
    ##############################################

    @contextmanager
    def component_start(self, name: str, description: str, skip: bool = False):
        """Start a new component block.
        The component will be automatically placed in the component tree.

        Args:
            name (str): Used to describe the component and for skipping by name
            description (str): Currently unused.
                Will be used to be a explanation of the overall goal of the component's checks
            skip (bool, optional): Controls whether all checks under the component should be skipped.
                This includes all child component checks. Defaults to False.
        """
        try:
            # set on protocol init
            # set current stage as parent
            new_component = ValidationProtocol._Component(
                name=name, parent=self.cur_component, description=description
            )
            self.cur_component = new_component

            directly_skipped_by_protocol_init = (
                new_component.name in self.skip_components
            )
            directly_skipped_by_add_call = skip
            parent_skipping_children = new_component.parent.skip_children
            not_in_run = not new_component.ancestry_is_in(self.run_components)
            explicitly_run_by_protocol_init = (
                new_component.name in self.run_components
                and self.run_components != ValidationProtocol._ALL_COMPONENTS()
            )

            # determine if skipping needs to be set for the new component
            if explicitly_run_by_protocol_init and not directly_skipped_by_add_call:
                pass
            elif any(
                [
                    directly_skipped_by_protocol_init,
                    directly_skipped_by_add_call,
                    parent_skipping_children,
                    not_in_run,
                ]
            ):
                new_component.skip = True
                # propagate skip children
                if any(
                    [
                        parent_skipping_children,
                        directly_skipped_by_protocol_init,
                        directly_skipped_by_add_call,
                    ]
                ):
                    new_component.skip_children = True
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
        """Provides a context to queue multiple checks with a common list of arguments.
        E.g. can be used to supply multiple files to a single 'check_file_exists' function.

        Args:
            payloads (list[dict]): A list of keyword arguments to pass to every check in the context

        """
        try:
            self._payloads = payloads
            yield
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
        """Adds the check to the queue for each payload.
        Payload can be either supplied directly on the add invocation
        or in a wrapping 'ValidationProtocol.payload' block.

        Args:
            fcn (Callable[..., FlagEntry]): The function to queue for each payload.
            payloads (dict, optional): A direct payload to supply to the fcn, defaults to None.
                Falls back to a wrapped 'ValidationProtocol.payload' if available
            skip (bool, optional): Denotes whether the queued check will be skipped. Defaults to False
            config (dict, optional): Additional function keyword arguments. Defaults to None.
                These are considered independent of the items that need validation (which should be supplied via payload)
            description (str, optional): A description of the check function. Defaults to function name.
                Should be used if the function name doesn't adequately describe what is being checked.
        """
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

    ##################################################################
    ### METHODS FOR DESCRIBING PLANNED VALIDATION CHECKS
    ##################################################################

    def queued_checks(
        self,
        include_individual_checks: bool = True,
        include_skipped_components: bool = False,
    ) -> str:
        """Returns a print-friendly string describing the queued checks.

        Args:
            include_individual_checks (bool, optional): Controls whether individual checks should be included. Defaults to True.

        Returns:
            str: A human friendly description of the queued checks.
        """
        INDENT_CHAR = " "

        # create by component dictionary
        check_by_component: dict[
            ValidationProtocol._Component, list[ValidationProtocol._QueuedCheck]
        ] = defaultdict(list)
        for check in self._check_queue:
            check_by_component[check["component"]].append(check)

        def sum_all_children(component: ValidationProtocol._Component) -> int:
            sum = len(check_by_component[component])
            for child in component.children:
                sum += sum_all_children(child)
            return sum

        # print tree of components from top down
        def render_self_and_children(component: ValidationProtocol._Component) -> str:
            if not include_skipped_components and component.skip:
                return ""
            INDENT_STR = INDENT_CHAR * (len(component.ancestor_line) - 1)
            count_str = f"[{sum_all_children(component)}"
            count_str2 = f"[{len(check_by_component[component])}"
            lead_str = f"{INDENT_STR}↳'{component.name}'{'-> !SKIPPED!' if component.skip else ''}"
            buffer = f"{lead_str : <55}DIRECT:{count_str2 : >4}] ALL:{count_str : >5}]"

            if include_individual_checks:
                check_lines = Counter(
                    [
                        f"{INDENT_STR} > {check['description']}"
                        for check in check_by_component[component]
                        if check["to_run"]
                    ]
                )
                check_line_print = [
                    f"{line} x {line_count}" for line, line_count in check_lines.items()
                ]
                if check_lines:
                    buffer += "\n" + "\n".join(check_line_print)

            for child in component.children:
                if render_self_and_children(child):
                    buffer += "\n" + render_self_and_children(child)
            return buffer

        return render_self_and_children(self._root_component)

    ##################################################################
    ### METHODS FOR RUNNING VALIDATION CHECKS
    ##################################################################
    def run(self, flag_unhandled_exceptions: bool = False):
        """Runs all queued checks not marked as skipped.

        Args:
            flag_unhandled_exceptions (bool, optional): If true, any unhandled exceptions within
                check functions will be flagged with a DEV_UNHANDLED flag but
                will not halt the generation of the overall protocol report.
                Defaults to False which means the protocol will raise the exception and stop midway.
        """
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
                try:
                    payload = self._eval_payload_callables(queued["payload"])
                except Exception as e:
                    raise RuntimeError(
                        f"Failed to evaluate payload: component:{queued['component']}, "
                        f"check function:{fcn_name}"
                    ) from e
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
                except Exception as e:
                    if flag_unhandled_exceptions:
                        packed_result = {
                            "code": FlagCode.DEV_UNHANDLED,
                            "message": str(e),
                            "function": fcn_name,
                            "description": queued["description"],
                            "kwargs": payload_and_config,
                            "config": queued["config"],
                        }
                    else:
                        raise RuntimeError(
                            f"Function failed: {fcn_name} part of {queued['component']}"
                        ) from e
            # add result (including skip flag) to component
            queued["component"].flags.append(packed_result)

    ##############################################
    ### METHODS FOR EXPORTING FLAG DATA
    ##############################################

    @staticmethod
    def _get_report_data(component):
        yield {tuple(component.ancestor_line): component.flags}
        for child in component.children:
            yield from ValidationProtocol._get_report_data(child)

    class Report(TypedDict):
        """A dictionary that acts as a container for sub-reports."""

        flag_table: pd.DataFrame
        """ A dataframe with each flag as a row. """
        outliers: pd.DataFrame
        """ A dataframe with outlier entities as rows and associated components as columns. """

    def report(
        self,
        include_skipped: bool = True,
    ) -> "ValidationProtocol.Report":
        """Tabulates the results of the executed protocol.

        Args:
            include_skipped (bool, optional): Controls whether the flag table should include skipped flags. Defaults to True.

        Returns:
            ValidationProtocol.Report: A report of all check results
        """
        COL_ORDER = [
            "description",
            "function",
            "code",
            "message",
            "code_level",
            "kwargs",
            "config",
        ]

        self.results = dict(self.results)
        data = list(self._get_report_data(self._root_component))
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

        # creating derivative colun from code
        df["code_level"] = df["code"].apply(lambda x: x.value)
        # resort columns
        df = df[COL_ORDER]

        return {
            "flag_table": df,
            "outliers": pd.DataFrame.from_dict(
                default_to_regular(self.outliers), orient="index"
            ).T,
        }

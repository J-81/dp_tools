from dp_tools.components.components import ReadsComponent
from dp_tools.core.check_model import Check, Flag, FlagCode


def _validate_func(self: Check, sample) -> Flag:
    if sample.dataset.metadata.paired_end:
        expected_components = ["rawForwardReads", "rawReverseReads"]
    else:
        expected_components = ["rawReads"]

    missing_components = list()
    unexpected_components = list()
    for expected_component in expected_components:
        component = getattr(sample, expected_component, None)
        if component == None:
            unexpected_components.append(expected_component)
        if not isinstance(component, ReadsComponent):
            missing_components.append(expected_component)

    if unexpected_components:
        code = FlagCode.DEV_HANDLED

    if missing_components:
        code = FlagCode.HALT1
    else:
        code = FlagCode.GREEN
    return Flag(
        check=self, code=code, message_args={"missing_components": missing_components}
    )


RAWREADS_0001 = Check(
    id="RAWREADS_0001",
    description=(
        "Check that appropriate raw reads components exist. Also check that "
        "All datafiles associated with the components are present. "
        "For paired end studies, this means both rawForwardReads and rawReverseReads "
        "Are attached components. For single end studies, "
        "this means the rawReads component is attached. "
    ),
    flag_desc={
        FlagCode.GREEN: "All expected raw read files present",
        FlagCode.HALT1: "Missing expected components: {missing_components}",
        FlagCode.DEV_HANDLED: "Searched for component, but component was not expected by entity model: {unexpected_components}",
    },
    validate_func=_validate_func,
)

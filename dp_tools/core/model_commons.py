""" Utilities common to all models in the package """

from typing import List
import logging
log = logging.getLogger(__name__)

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
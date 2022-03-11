#########################################################################
# FLAG (A validation report message)
#########################################################################
from dataclasses import dataclass
import enum
from typing import List

from dp_tools.core.model_commons import strict_type_checks


class FlagCode(enum.Enum):
    HALTING = 90
    RED = 50
    YELLOW = 30
    GREEN = 20
    INFO = 10

@dataclass
class Check:

    flags: List['Flag']

@dataclass
class Flag:

    code: FlagCode
    message: str
    check: Check
    
    # must ensure all flags are correct before continuing
    def __post_init__(self):
        strict_type_checks(self)



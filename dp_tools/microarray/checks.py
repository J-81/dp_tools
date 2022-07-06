from enum import Flag
from dp_tools.core.check_model import FlagCode, FlagEntry
import pandas as pd

def check_factor_values_in_runsheet(runsheet_file: str) -> FlagEntry: 
    """Verifies inclusion of factor values in runsheet

    :param runsheet_file: Path to runsheet file
    :type runsheet_file: str
    :return: determination if factor value was located
    :rtype: FlagEntry
    """
    runsheet = pd.read_csv(runsheet_file)
    code = FlagCode.HALT
    message = f"Runsheet does not have a Factor Value column"
    for col_name in runsheet.columns:
        if "Factor Value" in col_name:
            code = FlagCode.GREEN
            message = f"Factor Value column detected in runsheet"
    return {"code": code, "message": message}


    

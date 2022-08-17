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

def check_html_file(html_file: str) -> FlagEntry:
    """Verifies if the rendered Rmarkdown file has a valid html extension

    :param html_file: Path to rendered Rmarkdown
    :type html_file: str
    :return: determination if the file has a valid extension
    :rtype: FlagEntry
    """
    if html_file[-4:] == "html":
        code = FlagCode.GREEN
        message = f"File is an html file"
    else:
        code = FlagCode.HALT
        message = f"File is not an html file"
    return {"code": code, "message": message}






    

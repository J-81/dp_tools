from enum import Flag
from dp_tools.core.check_model import FlagCode, FlagEntry
import pandas as pd
#should I import data all in one place like the car in the demo?

def check_factor_values_in_runsheet(runsheet_file) -> FlagEntry: 
    """ Checks for Factor Value by iterating over column names
    :param file: Runsheet File
    :type file: Path
    :return: flagcode and message
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


# def check_file_exists(file: Path) -> FlagEntry:
#     # check logic
#     if file.is_file():
#         code = FlagCode.GREEN
#         message = f"File exists: {file.name} "
#     else:
#         code = FlagCode.HALT
#         message = f"Missing file: {file.name} expected at {str(file)} "

#     return {"code": code, "message": message}

def test_check_factor_values_in_runsheet():
    res = check_factor_values_in_runsheet(<filename>)
    assert res["code"] == FlagCode.Green
    assert res["message"]

# def main():
#     check_factor_values_in_runsheet(r"C:\Users\linde\rmarkdown\affy_processing\data_runsheets\GLDS_6_runsheet.csv")
# main()

#Possible checks (took from bulkRNASeq)
#if file exists 
#check for non null entries
#value content
#   what should the values look like?

    

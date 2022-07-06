from dp_tools.microarray.checks import check_factor_values_in_runsheet

from dp_tools.core.check_model import FlagCode


def test_check_factor_values_in_runsheet():
    res = check_factor_values_in_runsheet(r"C:\Users\linde\rmarkdown\affy_processing\data_runsheets\GLDS_6_runsheet.csv")
    assert res["code"] == FlagCode.GREEN
    assert res["message"]
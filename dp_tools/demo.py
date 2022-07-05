from dp_tools.core.check_model import ValidationProtocol, FlagCode

#########################################
# DATA MODEL INIT
#########################################

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


#########################################
# Check function definitions
#########################################


def check_if_wheel_flat(wheel: dict):
    if not wheel["isFlat"]:
        code = FlagCode.GREEN
        message = "Wheel is ready to go!"
    else:
        code = FlagCode.HALT
        message = "This wheel is flat!"
    return {"code": code, "message": message}


def check_engine_okay(car: dict):
    if car["engineOkay"]:
        code = FlagCode.GREEN
        message = "Engine looks good"
    else:
        code = FlagCode.HALT
        message = "Engine needs work!"
    return {"code": code, "message": message}


def check_gas(cur_gas: float, minimum_gas: float):
    if cur_gas >= minimum_gas:
        code = FlagCode.GREEN
        message = f"Gas tank is at {cur_gas}. Which is above {minimum_gas}"
    else:
        code = FlagCode.HALT
        message = (
            f"Gas tank needs a fill up! Current: {cur_gas}. Minimum: {minimum_gas}"
        )
    return {"code": code, "message": message}


def check_charge(cur_charge: float, minimum_charge: float):
    if cur_charge >= minimum_charge:
        code = FlagCode.GREEN
        message = "Charge looks good"
    else:
        code = FlagCode.HALT
        message = "Needs a charge!"
    return {"code": code, "message": message}


#########################################
# Protocol definition
#########################################

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


#########################################
# Running Protocol
#########################################

# Now printing the queued checks
print(vp.queued_checks())

print("********************************")
# Running the checks
vp.run()
vp.run(skip_components = ["Tires"])

# And showing the results
print(vp.report()["flag_table"])

print("Individual Flag Row")
print(vp.report()["flag_table"].iloc[5].to_dict())


#########################################
# Exporting Protocol Results to File
#########################################

vp.report()["flag_table"].to_csv("VV_log.tsv", sep = "\t")
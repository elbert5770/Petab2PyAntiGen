# Update_parameters.py: called once per replicate, after the model is built and before it is simulated, when the base
# parameters differ between replicates.  Weber's two conditions differ only in what happens on a clock (Events.py), so there is
# nothing to do: every value is in antimony_models/Weber_BMC2015/Weber_BMC2015_parameters.csv.

def update_no_parameters(r, replicate):
    return

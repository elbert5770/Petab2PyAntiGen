# Update_parameters.py: called once per replicate, after the model is built and before it is simulated.
# The conditions differ in five quantities, the columns of data/Smith_BMCSystBiol2013_conditions.csv that Experiment.py puts in
# the replicate. All five are parameters of the model (Ins and extracellular_ROS are boundary species in the SBML, held fixed;
# E2F1 is a species no reaction changes), so setting them is enough.

CONDITION_COLUMNS = ("Ins", "extracellular_ROS", "indicator_jnk", "indicator_foxo", "E2F1")


def update_Smith(r, replicate):
    for name in CONDITION_COLUMNS:
        r[name] = replicate[name]

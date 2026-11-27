# Update_parameters.py: called once per replicate, after the model is built and before it is simulated.
# The 19 conditions differ in the insulin they start with (data/Schwen_PONE2014_conditions.csv, put in the replicate by Experiment.py).
# It is set on the species, not on a parameter: PEtab's init_Ins is read by an initial assignment, which RoadRunner evaluates when
# the model is loaded.

def update_Schwen(r, replicate):
    r["Ins_ExtracellularMedium"] = replicate["init_Ins"]

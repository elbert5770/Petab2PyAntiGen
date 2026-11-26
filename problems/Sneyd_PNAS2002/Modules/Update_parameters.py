# Update_parameters.py: called once per replicate, after the model is built and before it is simulated.
# The conditions differ only in the (constant) IP3 and Ca2+, the two columns of data/Sneyd_PNAS2002_conditions.csv that
# Experiment.py puts in the replicate. Nothing in the model is an initial assignment of them, so setting the parameters is enough.

def update_Sneyd(r, replicate):
    r["IP3"] = replicate["IP3"]
    r["Ca"] = replicate["Ca"]

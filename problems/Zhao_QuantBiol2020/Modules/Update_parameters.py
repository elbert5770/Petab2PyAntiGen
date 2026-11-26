# Update_parameters.py: called once per replicate, after the model is built and before it is simulated.
# Each (region, stage) is its own replicate, and these are the values that distinguish them: the population, the initial
# infections and the rate parameters (PEtab's nominal ones; a fit replaces R, gamma_1 and gamma_2, see
# Update_opt_parameters.py).  They are the columns of data/Zhao_QuantBiol2020_conditions.csv, which Experiment.py puts in
# the replicate.
#
# The initial amounts are set on the species themselves, not through a parameter: RoadRunner evaluates a species'
# initial assignment when the model is loaded, so a parameter changed here would not reach it.

RATES = ("Total_Pop", "R", "gamma_1", "gamma_2")
INITIAL = {"Total_Pop": "Susceptible_Population", "U0": "Unquarantined_Infected_Population",
           "Q0": "Quarantined_Infected_Population", "C0": "Confirmed_Infected_Population"}


def update_Zhao(r, replicate):
    for name in RATES:
        r[name] = replicate[name]
    for column, species in INITIAL.items():
        r[species] = replicate[column]

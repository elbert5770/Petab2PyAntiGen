"""COVID-19 in Wuhan, Hubei and China, three stages of control (Zhao et al., Quantitative Biology 2020).

People move  Susceptible -> Unquarantined_Infected -> Quarantined_Infected -> Confirmed_Infected.
Unquarantined infected people infect the susceptible at rate  alpha * S * U / N  with  alpha = R * gamma_1; they are
quarantined at rate gamma_1; quarantined people are confirmed at rate beta = gamma_2 + (1 - gamma_2) * sigma.

PEtab selects the parameters, population and initial infections of each (region, stage) with six 0/1 "trigger"
parameters multiplied into the SBML rules, so every rate is a sum of seven products.  Here the model has the plain
parameters (R, gamma_1, gamma_2, Total_Pop, U0, Q0, C0) and each (region, stage) replicate sets its own values
(Modules/Experiment.py, Update_parameters.py).
"""
from pyantigen.generate.module_base import PyAntiGenModule

SPECIES = ["Susceptible", "Unquarantined_Infected", "Quarantined_Infected", "Confirmed_Infected"]


class SIQC_Zhao_QuantBiol2020(PyAntiGenModule):
    """Three reactions, and the rules that turn R, gamma_1, gamma_2 into the rates."""

    def build(self):
        Comp = self.config.get('Compartment', 'Population')
        S, U, Q, C = (f"{name}_{Comp}" for name in SPECIES)

        # infection: S -> U, catalysed by U; rate alpha * S * U / N
        self.add_reaction(f"Infection_{Comp}", f"[{S}]", f"[{U}]", "custom_conc_per_time",
                          f"alpha * {S} * {U} / Total_Pop")
        self.add_reaction(f"Quarantine_{Comp}", f"[{U}]", f"[{Q}]", "MA", "gamma_1")
        self.add_reaction(f"Confirmation_{Comp}", f"[{Q}]", f"[{C}]", "MA", "beta")

        self.add_rule("alpha := R * gamma_1")
        self.add_rule("beta := gamma_2 + (1 - gamma_2) * sigma")
        self.add_rule(f"Cumulative_Infected := {U} + {Q} + {C}")

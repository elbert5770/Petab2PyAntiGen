"""The type-2 IP3 receptor channel (Sneyd and Dufour, PNAS 2002), fitted to calcium release from liver microsomes.

Six states of one receptor, with first-order transitions whose rate constants depend on the calcium and IP3 of the
experiment (both constant within a condition):

    R  <-> O    activation by IP3 (and Ca2+)
    R  <-> I1   Ca2+ inhibition
    O  <-> A    the active state
    A  <-> I2   Ca2+ inhibition of the active state
    O  <-> S    the shut state

``L1, L3, L5`` are the dissociation constants that close the thermodynamic cycles; they are rules of the rate constants, not
parameters. The observed open probability is (0.9 A + 0.1 O)^4: four subunits, each open with that probability.
"""
from pyantigen.generate.module_base import PyAntiGenModule

STATES = ["IPR_O", "IPR_R", "IPR_I1", "IPR_S", "IPR_A", "IPR_I2"]
RATE_CONSTANTS = ["k1", "k2", "k3", "k4", "k_1", "k_2", "k_3", "k_4", "l2", "l4", "l6", "l_2", "l_4", "l_6"]

# (from, to, the first-order rate coefficient): each reaction is coefficient * [from]
TRANSITIONS = [
    ("IPR_O", "IPR_R",  "(k_2 + l_4 * Ca) / (1 + Ca / L5)"),
    ("IPR_R", "IPR_O",  "(k2 * L3 + l4 * Ca) / (L3 + Ca * (1 + L3 / L1)) * IP3"),
    ("IPR_R", "IPR_I1", "(k1 * L1 + l2) * Ca / (L1 + Ca * (1 + L1 / L3))"),
    ("IPR_I1", "IPR_R", "k_1 + l_2"),
    ("IPR_O", "IPR_A",  "(k4 * L5 + l6) * Ca / (L5 + Ca)"),
    ("IPR_A", "IPR_O",  "L1 * (k_4 + l_6) / (L1 + Ca)"),
    ("IPR_A", "IPR_I2", "(k1 * L1 + l2) * Ca / (L1 + Ca)"),
    ("IPR_I2", "IPR_A", "k_1 + l_2"),
    ("IPR_O", "IPR_S",  "k3 * L5 / (L5 + Ca)"),
    ("IPR_S", "IPR_O",  "k_3"),
]


class IP3Receptor_Sneyd_PNAS2002(PyAntiGenModule):
    """Ten transitions between six states, and the rules that close the thermodynamic cycles."""

    def build(self):
        Comp = self.config.get('Compartment', 'Membrane')
        for src, dst, coefficient in TRANSITIONS:
            self.add_reaction(f"{src}_to_{dst}_{Comp}", f"[{src}_{Comp}]", f"[{dst}_{Comp}]", "MA", f"({coefficient})")
        self.add_rule("L1 := k_1 / k1 / (l_2 / l2)")
        self.add_rule("L3 := k_2 / k2 / (l_4 / l4)")
        self.add_rule("L5 := k_4 / k4 / (l_6 / l6)")


class Observables_Sneyd_PNAS2002(PyAntiGenModule):
    """PEtab's one observable: the open probability."""

    def build(self):
        Comp = self.config.get('Compartment', 'Membrane')
        self.add_rule(f"open_probability := (0.9 * IPR_A_{Comp} + 0.1 * IPR_O_{Comp})^4")

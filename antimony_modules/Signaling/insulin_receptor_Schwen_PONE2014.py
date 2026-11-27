"""Insulin binding, internalisation and degradation in isolated hepatocytes (Schwen et al., PLoS One 2014).

Two subpopulations of cells, with low (1) and high (2) insulin binding, share one reaction pattern; the high-binding receptors differ by
the fold factors ``ka2fold`` and ``kd2fold`` on the association and dissociation constants (and have their own ``kin2``, ``kout2``).
For each receptor type ``i`` in (1, 2):

    Ins + Rec<i>  <->  IR<i>               binding (ka) and dissociation (kd)
    IR<i>          ->  IR<i>in             internalisation (kin)
    IR<i>in        ->  IR<i>               recycling (kout)
    IR<i>in        ->  Rec<i> + InsulinFragments   degradation (kout_frag): the receptor returns, the insulin leaves as fragments

and insulin also sticks to something else, ``Ins <-> BoundUnspec`` (kon_unspec, koff_unspec). The SBML's two ``Uptake`` species
(a running total, observed by nothing) and their reactions are dropped. Each experiment sets the insulin it starts with
(``init_Ins``, data/Schwen_PONE2014_conditions.csv); the receptor numbers are the initial amounts ``Rec1 = ini_R1`` and
``Rec2 = ini_R1 * ini_R2fold`` (InitialConditions.csv), so they are fitted through those two parameters.

PEtab's observables are measured on a log10 scale (``observableTransformation = log10``, normal noise): the rules here are the log10 of
a fitted scale times a signal, and the data file holds the log10 of the measurements. Insulin (ELISA) was measured in four experiment
series, each with its own offset, scale and saturation constant (``offset_nExpID<k>``, ``scaleElisa_nExpID<k>``, ``km_nExpID<k>``): one rule
per series.
"""
from pyantigen.generate.module_base import PyAntiGenModule

COMP = "ExtracellularMedium"
EXPERIMENTS = (1, 2, 3, 4)          # the four ELISA series

# receptor type -> (binding, dissociation, internalisation, recycling) rate constants, as a rate expression
RATES = {1: ("ka1", "kd1", "kin", "kout"),
         2: ("ka1 * ka2fold", "kd1 * kd2fold", "kin2", "kout2")}


def sp(name):
    return f"{name}_{COMP}"


class InsulinReceptor_Schwen_PONE2014(PyAntiGenModule):
    """The 12 reactions: five per receptor type and the unspecific binding."""

    def build(self):
        Ins, Frag, Bound = sp("Ins"), sp("InsulinFragments"), sp("BoundUnspec")

        def rxn(name, reactants, products, rate):
            self.add_reaction(name, " + ".join(f"[{s}]" for s in reactants), " + ".join(f"[{s}]" for s in products), "MA", rate)

        for i, (ka, kd, kin, kout) in RATES.items():
            Rec, IR, IRin = sp(f"Rec{i}"), sp(f"IR{i}"), sp(f"IR{i}in")
            rxn(f"Binding{i}", [Ins, Rec], [IR], ka)
            rxn(f"Dissociation{i}", [IR], [Ins, Rec], kd)
            rxn(f"Internalisation{i}", [IR], [IRin], kin)
            rxn(f"Recycling{i}", [IRin], [IR], kout)
            rxn(f"Degradation{i}", [IRin], [Rec, Frag], "kout_frag")
        rxn("UnspecificBinding", [Ins], [Bound], "kon_unspec")
        rxn("UnspecificRelease", [Bound], [Ins], "koff_unspec")


class Observables_Schwen_PONE2014(PyAntiGenModule):
    """PEtab's four observables as log10 rules: IR1, IR2 and IRsum (cytometry) and insulin (ELISA, one rule per series)."""

    def build(self):
        IR1, IR2, IR1in, IR2in = (sp(s) for s in ("IR1", "IR2", "IR1in", "IR2in"))
        Ins, Frag = sp("Ins"), sp("InsulinFragments")
        self.add_rule(f"observable_IR1 := log10(scale * ({IR1} + {IR1in} + offset))")
        self.add_rule(f"observable_IR2 := log10(scale * ({IR2} + {IR2in} + offset))")
        self.add_rule(f"observable_IRsum := log10(scale * (0.605 * {IR1} + 0.395 * {IR2} + 0.605 * {IR1in} + 0.395 * {IR2in} + offset))")
        for k in EXPERIMENTS:
            total = f"({Ins} + {Frag} * fragments)"            # the ELISA also sees insulin fragments, to the extent ``fragments``
            self.add_rule(f"observable_Insulin_nExpID{k} := log10(offset_nExpID{k} + scaleElisa_nExpID{k} * {total} / ({total} / km_nExpID{k} + 1))")

"""Lipid-transfer signalling at the Golgi: PKD, PI4KIIIbeta and CERT (Weber et al., BMC Systems Biology 2015).

Three kinases/transfer proteins, each cycling between an inactive and an active form, with feedback between them:

    PKD   <-> PKDDAGa     activated by DAG (PdBu mimics it, u5), by CERTERa/PI4K3Ba, and inactivated (inhibitor, u6)
    PI4K3B <-> PI4K3Ba    activated by PKDDAGa
    CERT -> CERTERa -> CERTTGNa -> CERT   shuttling ER -> trans-Golgi, driven by PI4K3Ba and PKDDAGa

Every species is made at a constant rate and degraded first-order. Four time-dependent inputs u3 (ectopic PI4K3B), u4 (ectopic
CERT), u5 (PdBu) and u6 (the PKD inhibitor NB142-70) are 0 until an event sets them (Modules/Events.py): PEtab wrote each as a
SBML rule ``piecewise(0, time < t_on, dose)`` driven by condition parameters. PEtab's observables, a signal times a fitted scale
factor, are rules here (Observables_Weber_BMC2015); the scale factors are ordinary parameters.
"""
from pyantigen.generate.module_base import PyAntiGenModule

# species the module creates, for the observables and the data
SPECIES = ["PKD", "PKDDAGa", "PI4K3B", "PI4K3Ba", "CERTERa", "CERT", "CERTTGNa"]


class SignalingCascade_Weber_BMC2015(PyAntiGenModule):
    """The 20 reactions. A rate is a concentration per hour; the compartment volume V_Cyt is applied by the generator."""

    def build(self):
        Comp = self.config.get('Compartment', 'Cyt')
        PKD, PKDa, PI, PIa, CE, C, CT = (f"{s}_{Comp}" for s in SPECIES)

        def rxn(name, reactant, product, rate_type, rate):
            self.add_reaction(f"{name}_{Comp}", f"[{reactant}]" if reactant else "", f"[{product}]" if product else "",
                              rate_type, rate)

        # PKD <-> PKDDAGa
        rxn("PKD_activation_by_CERT_PI4K3B", PKD, PKDa, "custom_conc_per_time",
            f"{CE} * {PIa} * {PKD} * p11 * p31 / (({PIa} + m31) * (m11 + {CE} * {PIa} * p31 / ({PIa} + m31)))")
        rxn("PKD_activation_by_DAG", PKD, PKDa, "custom_conc_per_time", f"{PKD} * p12 * (pu5 * u5 + 1)")
        rxn("PKD_inactivation", PKDa, PKD, "custom_conc_per_time", f"{PKDa} * p13 * (pu6 * u6 + 1)")
        # PI4K3B <-> PI4K3Ba
        rxn("PI4K3B_inactivation", PIa, PI, "MA", "p21")
        rxn("PI4K3B_activation_by_PKD", PI, PIa, "custom_conc_per_time", f"{PI} * {PKDa} * p22 / ({PKDa} + m22)")
        # CERT -> CERTERa -> CERTTGNa -> CERT
        rxn("CERT_to_ER", C, CE, "MA", "p32")
        rxn("CERT_ER_to_TGN", CE, CT, "custom_conc_per_time", f"{CE} * {PIa} * p31 / ({PIa} + m31)")
        rxn("CERT_TGN_release", CT, C, "custom_conc_per_time", f"{CT} * {PKDa} * p33 / ({PKDa} + m33)")

        # constant synthesis; the ectopic expression of PI4K3B and CERT is the input u3 / u4
        rxn("PKD_synthesis", "", PKD, "MA", "s12")
        rxn("PI4K3B_synthesis", "", PI, "MA", "s21")
        rxn("PI4K3B_ectopic_expression", "", PI, "MA", "pu3 * u3")
        rxn("CERTERa_synthesis", "", CE, "MA", "s31")
        rxn("CERTERa_ectopic_expression", "", CE, "MA", "pu4 * u4")

        # first-order degradation
        for species, rate in zip((PKD, PKDa, PI, PIa, CE, C, CT), ("a11", "a12", "a21", "a22", "a31", "a32", "a33")):
            rxn(f"Degradation_{species.split('_')[0]}", species, "", "MA", rate)


class Observables_Weber_BMC2015(PyAntiGenModule):
    """PEtab's eight observables: three totals and five signals, each the active fraction (or amount) times a fitted scale."""

    def build(self):
        Comp = self.config.get('Compartment', 'Cyt')
        PKD, PKDa, PI, PIa, CE, C, CT = (f"{s}_{Comp}" for s in SPECIES)
        self.add_rule(f"yPKDt := {PKD} + {PKDa}")
        self.add_rule(f"yPI4K3Bt := {PI} + {PIa}")
        self.add_rule(f"yCERTt := {C} + {CE} + {CT}")
        # phospho-signals: scale * active amount (PKD), scale * active fraction (PI4K3B, CERT)
        for name in ("yPKDpN0", "yPKDpN24", "yPKDpN25"):
            self.add_rule(f"{name} := {PKDa} * scale_{name}")
        self.add_rule(f"yPI4K3BpRN24 := {PIa} * scale_yPI4K3BpRN24 / ({PI} + {PIa})")
        self.add_rule(f"yCERTpRN24 := {C} * scale_yCERTpRN24 / ({C} + {CE} + {CT})")

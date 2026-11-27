"""Insulin signalling, oxidative stress and FOXO1 in an adipocyte (Smith et al., BMC Systems Biology 2013).

Four blocks, all with amounts as species and rates in amount/time (rate type ``custom``: the expression is used as written, species are amounts):

    InsulinSignaling_Smith_BMCSystBiol2013   insulin receptor (binding, phosphorylation, internalisation), IRS1, PI3K, PIP2/PIP3,
                                              Akt, PKC, AS160 and the GLUT4 transporter that carries glucose uptake
    OxidativeStress_Smith_BMCSystBiol2013    ROS made by NOX and the mitochondria, removed by SOD2 and GSH; ROS oxidises PTP1B, PTEN
                                              and DUSP; JNK and IKK are switched on by ROS (and off by DUSP)
    Foxo1_Smith_BMCSystBiol2013               FOXO1 in 16 modification states (Pa, Pd, Pe: phosphorylated by Akt/SGK, IKK, JNK; pUb:
                                              ubiquitinated) in three places (cytoplasm, nucleus, DNA-bound): 48 species from one loop;
                                              transcription of InR and SOD2, their mRNA and protein
    Observables_Smith_BMCSystBiol2013        PEtab's nine observables, each a fitted scale times a model quantity

The FOXO1 rates depend on the modification bits through a few named factors (``f_Pe_cn``: the factor by which phosphorylation Pe
changes the cytoplasm -> nucleus rate), instead of PEtab's 150 separate rate constants. ``Ins``, ``extracellular_ROS`` and ``E2F1``
are parameters that each experiment sets (data/Smith_BMCSystBiol2013_conditions.csv); ``Ins`` and ``extracellular_ROS`` are
boundary species in the SBML, held fixed. ``indicator_jnk`` and ``indicator_foxo`` switch the JNK/IKK and FOXO1 blocks off in the
insulin-dose experiments.
"""
import itertools
import re

from pyantigen.generate.module_base import PyAntiGenModule

CYT, CS, NUC, DNA = "Cytoplasm", "Cellsurface", "Nucleus", "Dnabound"
MODS = ("Pa", "Pd", "Pe", "pUb")                  # the four FOXO1 modifications, each 0 or 1
STATES = list(itertools.product((0, 1), repeat=4))


_SPECIES = re.compile(r"\b\w+?_(Cytoplasm|Cellsurface|Nucleus|Dnabound)\b")


def amounts(rate):
    """PyAntiGen reads a species in a rate law as a concentration, ``X`` -> ``X / V_<Compartment>``. The laws of this model are in
    amounts (as in the SBML, whose compartment sizes only enter once, in the ROS permeation), so each species is written as
    ``V_<Compartment> * X``: the generated rate is the SBML's, and the module reads like it."""
    return _SPECIES.sub(lambda m: m.group(0) if m.group(0).startswith("V_") else f"(V_{m.group(1)} * {m.group(0)})", rate)


def foxo(state, place):
    """Name of the FOXO1 species with modification bits ``state`` (Pa, Pd, Pe, pUb) in ``place``."""
    return "Foxo1_" + "_".join(f"{m}{b}" for m, b in zip(MODS, state)) + f"_{place}"


def foxo_species():
    return [foxo(s, place) for place in (CYT, NUC, DNA) for s in STATES]


def _factors(state, route, mods):
    """``f_<mod>_<route>`` for every modification in ``mods`` that is on in ``state``: the product scales a base rate."""
    return "".join(f" * f_{m}_{route}" for m in mods if state[MODS.index(m)])


class InsulinSignaling_Smith_BMCSystBiol2013(PyAntiGenModule):
    """Config: ``pip3_switch_eps`` (default None: PEtab's exact ``piecewise`` switch of Akt and PKC)."""

    def build(self):
        def rxn(name, reactants, products, rate):
            self.add_reaction(name, reactants, products, "custom", amounts(rate))

        c = lambda s: f"{s}_{CYT}"
        s_ = lambda s: f"{s}_{CS}"
        InR, InsInR, InsInRP, Ins2InRP = (s_(x) for x in ("InR", "Ins_InR", "Ins_InR_P", "Ins_2_InR_P"))
        InR_c, Ins2InRP_c, InsInRP_c = (c(x) for x in ("InR", "Ins_2_InR_P", "Ins_InR_P"))
        PTP1B, PP2A = c("PTP1B"), "PP2A"
        IRS1, IRS1P, IRS1S, IRS1PS = (c(x) for x in ("IRS1", "IRS1_TyrP", "IRS1_PolySerP", "IRS1_TyrP_PolySerP"))

        # receptor: insulin binds twice, the receptor is phosphorylated, dephosphorylated by PTP1B, and cycles through the cell
        rxn("R1f", f"[{InR}]", f"[{InsInR}]", f"k1 * Ins * {InR}")
        rxn("R1r", f"[{InsInR}]", f"[{InR}]", f"kminus1 * {InsInR}")
        rxn("R2", f"[{InsInR}]", f"[{InsInRP}]", f"k3 * {InsInR}")
        rxn("R3f", f"[{InsInRP}]", f"[{Ins2InRP}]", f"k2 * Ins * {InsInRP}")
        rxn("R3r", f"[{Ins2InRP}]", f"[{InsInRP}]", f"kminus2 * {Ins2InRP}")
        rxn("R4", f"[{InsInRP}]", f"[{InR}]", f"kminus3 * {PTP1B} * {InsInRP}")
        rxn("R5f", f"[{InR}]", f"[{InR_c}]", f"k4 * {InR}")
        rxn("R5r", f"[{InR_c}]", f"[{InR}]", f"kminus4 * {InR_c}")
        for name, surface, inside in (("R6", Ins2InRP, Ins2InRP_c), ("R7", InsInRP, InsInRP_c)):
            rxn(f"{name}f", f"[{surface}]", f"[{inside}]", f"k4prime * {surface}")
            rxn(f"{name}r", f"[{inside}]", f"[{surface}]", f"kminus4prime * {inside}")
        for name, inside in (("R10", Ins2InRP_c), ("R11", InsInRP_c)):
            rxn(name, f"[{inside}]", f"[{InR_c}]", f"k6 * {PTP1B} * {inside}")

        # IRS1: made, degraded, phosphorylated on tyrosine by the receptor (and back), on serine by PKC, IKK and JNK; binds PI3K
        PKC_P, IKK_P, JNK_P = (c(x) for x in ("PKC_P", "IKK_P", "JNK_P"))
        PI3K, IRS1PI3K = c("PI3K"), c("IRS1_TyrP_PI3K")
        rxn("R12f", f"[{IRS1}]", f"[{IRS1P}]", f"cyto_vol * (k7 * {IRS1} * ({Ins2InRP} + {InsInRP}) / IRp)")
        rxn("R12r", f"[{IRS1P}]", f"[{IRS1}]", f"cyto_vol * kminus7a * {PTP1B} * {IRS1P}")
        for name, free, bound in (("R12_a", IRS1, IRS1S), ("R12_b", IRS1P, IRS1PS)):
            rxn(f"{name}_f", f"[{free}]", f"[{bound}]", f"cyto_vol * k2psp * {free} * {PKC_P}")
            rxn(f"{name}_r", f"[{bound}]", f"[{free}]", f"cyto_vol * kminus7b * {PP2A} * {bound}")
        rxn("R13f", f"[{PI3K}] + [{IRS1P}]", f"[{IRS1PI3K}]", f"cyto_vol * k8 * {IRS1P} * {PI3K}")
        rxn("R13r", f"[{IRS1PI3K}]", f"[{PI3K}] + [{IRS1P}]", f"cyto_vol * kminus8 * {IRS1PI3K}")
        rxn("R50f", "", f"[{IRS1}]", "k_irs1_basal_syn")
        for name, species, k in (("R50r1", IRS1, "k_irs1_basal_degr"), ("R50r2", IRS1P, "k_irs1_basal_degr"),
                                 ("R50r3", IRS1S, "k_irs1_polyserp_degr"), ("R50r4", IRS1PS, "k_irs1_polyserp_degr")):
            rxn(name, f"[{species}]", "", f"{k} * {species}")
        rxn("R51f", f"[{IRS1}]", f"[{IRS1S}]", f"kcat51 * {IRS1} * {IKK_P} / (Km51 + {IRS1})")
        rxn("R52f", f"[{IRS1}]", f"[{IRS1S}]", f"kcat52 * {IRS1} * {JNK_P} / (Km52 + {IRS1})")

        # PIP2 <-> PIP3 (made by IRS1-bound PI3K, removed by PTEN); Akt and PKC switch on above a PIP3 threshold; PP2A switches off
        PIP2, PIP3, PTEN = c("PIP2"), c("PI345P3"), c("PTEN")
        rxn("R14f", f"[{PIP2}]", f"[{PIP3}]", f"(k9_basal + k9 * {IRS1PI3K}) * {PIP2}")
        rxn("R14r", f"[{PIP3}]", f"[{PIP2}]", f"(kminus9_basal + kminus9 * {PTEN}) * {PIP3}")
        # PEtab: piecewise(x, x > 0, 0) = max(x, 0), with x = PIP3 - threshold. AMICI integrates a piecewise as an event, and where PIP3
        # crosses the threshold CVODE sometimes fails (see Optimizer_settings.py). ``pip3_switch_eps`` (molecules) replaces it by the
        # smooth (x + sqrt(x^2 + eps^2)) / 2, which differs from max(x, 0) by at most eps / 2: no event, but no longer PEtab's model.
        eps = self.config.get("pip3_switch_eps")
        x = f"({PIP3} - pip3_basal)"
        above = f"piecewise({x}, {PIP3} > pip3_basal, 0)" if not eps else f"(({x} + sqrt({x}^2 + {eps}^2)) / 2)"
        for name, off, on, kf, kr in (("R16", c("Akt"), c("Akt_P2"), "k11", "kminus11"), ("R17", c("PKC"), PKC_P, "k12", "kminus12")):
            rxn(f"{name}f", f"[{off}]", f"[{on}]", f"{kf} * {off} * {above}")
            rxn(f"{name}r", f"[{on}]", f"[{off}]", f"{kr} * {PP2A} * {on}")
        AS160, AS160_P = c("AS160"), c("AS160_P")
        rxn("R16a_f", f"[{AS160}]", f"[{AS160_P}]", f"kr16a * {c('Akt_P2')} * {AS160}")
        rxn("R16a_r", f"[{AS160_P}]", f"[{AS160}]", f"kminusr16a * {PP2A} * {AS160_P}")
        # GLUT4 to the cell surface, driven by phosphorylated AS160: what glucose uptake measures
        rxn("R20f", f"[{c('GLUT4')}]", f"[{s_('GLUT4')}]", f"(k13_basal + k13 * {AS160_P}) * {c('GLUT4')}")
        rxn("R20r", f"[{s_('GLUT4')}]", f"[{c('GLUT4')}]", f"kminus13 * {s_('GLUT4')}")


class OxidativeStress_Smith_BMCSystBiol2013(PyAntiGenModule):
    def build(self):
        def rxn(name, reactants, products, rate):
            self.add_reaction(name, reactants, products, "custom", amounts(rate))

        c = lambda s: f"{s}_{CYT}"
        ROS, GSH, GSSG, SOD2 = (c(x) for x in ("ROS", "GSH", "GSSG", "SOD2"))
        NOX_inact, NOX, NOX_deact = (c(x) for x in ("NOX_inact", "NOX", "NOX_deact"))

        # PTP1B and PTEN: oxidised by ROS, reduced by GSH.  (The enzymes are catalysts: ROS and GSH are not used up.)
        for name, enzyme, kf, kr in (("R30", "PTP1B", "k30f", "k30r"), ("R31", "PTEN", "k31f", "k31r")):
            rxn(f"{name}f", f"[{c(enzyme)}]", f"[{c(enzyme + '_ox')}]", f"{kf} * {c(enzyme)} * {ROS}")
            rxn(f"{name}r", f"[{c(enzyme + '_ox')}]", f"[{c(enzyme)}]", f"{kr} * {c(enzyme + '_ox')} * {GSH}")
        # NOX: activated by insulin, inactivates itself; makes ROS.  SOD2 removes ROS; mitochondria make it at a constant rate
        rxn("R34f", f"[{NOX_inact}]", f"[{NOX}]", f"k34f * {NOX_inact} * Ins")
        rxn("R34r1", f"[{NOX}]", f"[{NOX_deact}]", f"k34r1 * {NOX} * {NOX}")
        rxn("R34r2", f"[{NOX}]", f"[{NOX_inact}]", f"k34r2 * {NOX}")
        rxn("R34r3", f"[{NOX_deact}]", f"[{NOX_inact}]", f"k34r3 * {NOX_deact} * {NOX_deact}")
        rxn("R35f", "", f"[{ROS}]", f"k35f * {NOX}")
        rxn("R35r", f"[{ROS}]", "", f"k35r * {ROS} * {SOD2}")
        rxn("R36f", "", f"[{ROS}]", "k36f * Mt")
        # H2O2 crosses the membrane both ways, with the same permeability
        rxn("R37f", "", f"[{ROS}]", "ros_perm * membrane_area * extracellular_ROS")
        rxn("R37r", f"[{ROS}]", "", f"ros_perm * membrane_area * (V_Extracellular / V_Cytoplasm) * {ROS}")
        rxn("R38f", f"[{GSH}]", f"[{GSSG}]", f"k38f * {GSH} * {ROS}")
        rxn("R38r", f"[{GSSG}]", f"[{GSH}]", f"k38r * {GSSG}")

        # JNK and IKK are phosphorylated by ROS and dephosphorylated by DUSP, which ROS oxidises and GSH reduces
        DUSP, DUSP_ox = c("DUSP"), c("DUSP_ox")
        for name, kinase, kf, kr in (("R42", "JNK", "k42f * alpha_ox", "k42r"), ("R43", "IKK", "k43f", "k43r")):
            rxn(f"{name}f", f"[{c(kinase)}]", f"[{c(kinase + '_P')}]", f"{kf} * {c(kinase)} * {ROS} * indicator_jnk")
            rxn(f"{name}r", f"[{c(kinase + '_P')}]", f"[{c(kinase)}]", f"{kr} * {c(kinase + '_P')} * {DUSP} * indicator_jnk")
        rxn("R32f", f"[{DUSP}]", f"[{DUSP_ox}]", f"k32f * {DUSP} * {ROS} * indicator_jnk")
        rxn("R32r", f"[{DUSP_ox}]", f"[{DUSP}]", f"k32r * {DUSP_ox} * {GSH} * indicator_jnk")


class Foxo1_Smith_BMCSystBiol2013(PyAntiGenModule):
    """FOXO1: 16 modification states x (cytoplasm, nucleus, DNA-bound), and the two genes it drives (InR and SOD2)."""

    # (from, to, route): the four transport steps, and the modifications that change each one's rate
    TRANSPORT = ((CYT, NUC, "cn", ("Pe", "Pd", "Pa")), (NUC, CYT, "nc", ("Pe", "Pd", "Pa")),
                 (NUC, DNA, "nd", ("Pa",)), (DNA, NUC, "dn", ()))
    # modification -> (the kinase that adds it, its rate constant)
    KINASE = {"Pa": ("(Akt_P2_Cytoplasm + SGK)", "kkin_Pa"), "Pd": ("IKK_P_Cytoplasm * by_ikk_phos_factor", "kkin_Pd"),
              "Pe": ("JNK_P_Cytoplasm * by_jnk_phos_factor", "kkin_Pe")}
    GENES = ("InR", "SOD2")

    def build(self):
        def rxn(name, reactants, products, rate):
            self.add_reaction(name, reactants, products, "custom", amounts(rate))

        # synthesis of unmodified FOXO1, driven by E2F1
        rxn("Foxo_synthesis", "", f"[{foxo((0, 0, 0, 0), CYT)}]", "E2F1 * kfoxo_synth * indicator_foxo")
        for state in STATES:
            for src, dst, route, mods in self.TRANSPORT:
                rxn(f"Foxo_{route}_{''.join(map(str, state))}", f"[{foxo(state, src)}]", f"[{foxo(state, dst)}]",
                    f"{foxo(state, src)} * ktr_{route}{_factors(state, route, mods)} * indicator_foxo")
            for place in (CYT, NUC, DNA):
                for i, mod in enumerate(MODS[:3]):
                    on = tuple(1 if j == i else b for j, b in enumerate(state))
                    off = tuple(0 if j == i else b for j, b in enumerate(state))
                    if state[i] == 0:                                  # add the modification
                        kinase, k = self.KINASE[mod]
                        rxn(f"Foxo_{mod}_on_{place}_{''.join(map(str, state))}", f"[{foxo(state, place)}]", f"[{foxo(on, place)}]",
                            f"{foxo(state, place)} * {kinase} * {k} * indicator_foxo")
                    else:                                              # PP2A takes it off
                        rxn(f"Foxo_{mod}_off_{place}_{''.join(map(str, state))}", f"[{foxo(state, place)}]", f"[{foxo(off, place)}]",
                            f"{foxo(state, place)} * PP2A * kdephos * indicator_foxo")
                if state[3] == 0:                                      # ubiquitination by SCF (rate set by Pa and Pd)
                    rxn(f"Foxo_ub_{place}_{''.join(map(str, state))}", f"[{foxo(state, place)}]",
                        f"[{foxo(state[:3] + (1,), place)}]",
                        f"{foxo(state, place)} * SCF * kub{_factors(state, 'ub', ('Pd', 'Pa'))} * indicator_foxo")
            if state[3] == 1:                                          # ubiquitinated cytoplasmic FOXO1 is degraded
                rxn(f"Foxo_degradation_{''.join(map(str, state))}", f"[{foxo(state, CYT)}]", "",
                    f"{foxo(state, CYT)} * Proteasome * kdeg * indicator_foxo")

        # DNA-bound FOXO1, in any state, drives the transcription of InR and SOD2 (on top of a basal rate): mRNA is exported from the
        # nucleus, degraded, and translated into the protein, which is degraded
        bound = " + ".join(foxo(s, DNA) for s in STATES)
        for g in self.GENES:
            rna_nuc, rna_cyt, protein = f"RNA_{g}_{NUC}", f"RNA_{g}_{CYT}", f"{g}_{CYT}"
            rxn(f"{g}_transcription", "", f"[{rna_nuc}]", f"({bound}) * ktranscr_{g} * indicator_foxo")
            rxn(f"{g}_basal_transcription", "", f"[{rna_nuc}]", f"kbasal_{g} * indicator_foxo")
            rxn(f"{g}_mRNA_export", f"[{rna_nuc}]", f"[{rna_cyt}]", f"{rna_nuc} * kexp * indicator_foxo")
            rxn(f"{g}_mRNA_degradation", f"[{rna_cyt}]", "", f"{rna_cyt} * kmdeg * indicator_foxo")
            rxn(f"{g}_translation", "", f"[{protein}]", f"{rna_cyt} * ktransl_{g} * indicator_foxo")
            rxn(f"{g}_degradation", f"[{protein}]", "", f"{protein} * kpdeg_{g} * indicator_foxo")


class Observables_Smith_BMCSystBiol2013(PyAntiGenModule):
    """PEtab's nine observables: a fitted scale times a model quantity. All have noise sigma = 1."""

    def build(self):
        c = lambda s: f"{s}_{CYT}"
        glut4 = f"GLUT4_{CS}"
        bound_insulin = " + ".join(f"{x}_{CS}" for x in ("Ins_InR", "Ins_InR_P", "Ins_2_InR_P"))
        foxo_all = " + ".join(foxo_species())
        for name, formula in (
                ("PI3K_activity__2B", f"sc_PI3K * {c('IRS1_TyrP_PI3K')}"),
                ("Glucose_uptake__2B", f"sc_GLUT_2B * {glut4}"),
                ("Cell_Bound_Ins__2B", f"sc_BINS * ({bound_insulin})"),
                ("IRSYp__2C", f"sc_PIRS * {c('IRS1_TyrP')}"),
                ("PTP_activ__2D", f"sc_PTP * {c('PTP1B')} / ({c('PTP1B')} + {c('PTP1B_ox')})"),
                ("Glucose_uptake__120__3B", f"sc_GLUT_3B_120 * {glut4}"),
                ("Glucose_uptake__240__3B", f"sc_GLUT_3B_240 * {glut4}"),
                ("MnSOD_fold_induction__3C", f"sc_SOD2 * {c('SOD2')}"),
                ("FOXO4__3C", f"sc_FOXO1 * ({foxo_all})")):
            self.add_rule(f"{name} := {formula}")

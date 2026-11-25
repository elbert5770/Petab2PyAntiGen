"""Combinatorial methylation of histone H3 lysines 27 and 36 (Zheng et al., PNAS 2012).

An H3 tail carries 0-3 methyl groups on K27 and 0-3 on K36.  The 4 x 4 grid of
states K27me<a>K36me<b> is a lattice: one methyl group is added to, or removed
from, one lysine at a time.  K27me3K36me3 is not observed and is not in the
model, which leaves 15 states and 22 neighbouring pairs.  Newly made H3 enters
unmodified, and every state is diluted at the same rate by that inflow
(``inflow = inflowp * dilution``; the experiment switches dilution off).

PEtab spells this out as 60 reactions, each wrapped in its own SBML function
definition.  Here it is three loops.
"""
from pyantigen.generate.module_base import PyAntiGenModule

N_ME = 4                                                       # me0 .. me3 per lysine
STATES = [(a, b) for a in range(N_ME) for b in range(N_ME) if (a, b) != (3, 3)]


def state(a, b):
    """Species stem of K27me<a>K36me<b>."""
    return f"K27me{a}K36me{b}"


def steps():
    """Every methylation step as ((a, b), (c, d)): one methyl group added on K27 or K36."""
    return [((a, b), (c, d)) for (a, b) in STATES for (c, d) in [(a + 1, b), (a, b + 1)]
            if (c, d) in STATES]


def rate_constants():
    """Names of the 44 first-order rate constants, forward then reverse for each step."""
    return [k for (a, b), (c, d) in steps() for k in (f"k{a}{b}_{c}{d}", f"k{c}{d}_{a}{b}")]


class MethylationGrid_Zheng_PNAS2012(PyAntiGenModule):
    """Reversible methylation steps, inflow of unmodified H3, and dilution of every state."""

    def build(self):
        Comp = self.config.get('Compartment', 'Default')
        for Isotope in self.model.isotopes:
            Isotope_str = f"_{Isotope}_" if Isotope else "_"

            def sp(a, b):
                return f"{state(a, b)}{Isotope_str}{Comp}"

            # add one methyl group on K27 (a -> a+1) or on K36 (b -> b+1); one RMA
            # reaction carries both directions: Rate_eqtn_prototype = [k_forward, k_reverse]
            for (a, b), (c, d) in steps():
                self.add_reaction(
                    f"Methylation_{state(a, b)}_{state(c, d)}{Isotope_str}{Comp}",
                    f"[{sp(a, b)}]", f"[{sp(c, d)}]", "RMA",
                    f"[k{a}{b}_{c}{d}, k{c}{d}_{a}{b}]")

            # newly synthesised, unmodified H3 ...
            self.add_reaction(f"Inflow_{state(0, 0)}{Isotope_str}{Comp}",
                              "", f"[{sp(0, 0)}]", "MA", "inflow")
            # ... and dilution of every state at the same rate
            for (a, b) in STATES:
                self.add_reaction(f"Dilution_{state(a, b)}{Isotope_str}{Comp}",
                                  f"[{sp(a, b)}]", "", "MA", "inflow")

        self.add_rule("inflow := inflowp * dilution")


class Observables_Zheng_PNAS2012(PyAntiGenModule):
    """Each state as a fraction of all H3, the quantity PEtab's 15 observables measure."""

    def build(self):
        Comp = self.config.get('Compartment', 'Default')
        total = " + ".join(f"{state(a, b)}_{Comp}" for a, b in STATES)
        self.add_rule(f"H3_total := {total}")
        for (a, b) in STATES:
            self.add_rule(f"observable_{state(a, b)} := {state(a, b)}_{Comp} / H3_total")

"""How much does PyAntiGen compress Zheng_PNAS2012?  Two questions, measured separately.

  A. THE MODEL.  SBML -> Antimony (Tellurium) -> the Antimony PyAntiGen generates -> the module that generates it.
     Does the module say less than the Antimony it produces, and how does each grow with the size of the system?
  B. THE ESTIMATION PROBLEM.  PEtab's parameters / observables / conditions / measurements / problem tables -> PyAntiGen's
     Optimizer_settings / Loss_config / Events + Solver_settings / Data / Experiment.

    python benchmarks/compression_zheng.py [--petab <dir with problem.yaml>]

Writes benchmarks/COMPRESSION.md and benchmarks/compression.png.  Sizes are measured four ways (petab2pyantigen.compression):
bytes, non-blank lines, code lines (no comments or docstrings), and gzip -9, the baseline a representation has to beat to
claim it removed information and not only repetition.
"""
import argparse
import math
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from petab2pyantigen.compression import (Size, class_range, function_range, measure, measure_text,  # noqa: E402
                                         total)
from petab2pyantigen.petab_reader import load_petab                                             # noqa: E402
from petab2pyantigen.sbml_to_antimony import sbml_to_antimony                                   # noqa: E402

MODEL = "Zheng_PNAS2012"
MODULE = ROOT / "antimony_modules" / "Histone" / "methylation_Zheng_PNAS2012.py"
PROJECT = ROOT / "problems" / MODEL
MODS = PROJECT / "Modules"
AM = ROOT / "antimony_models" / MODEL


# --------------------------------------------------------------------------------------------------------------------
# generated Antimony, for any lattice size
# --------------------------------------------------------------------------------------------------------------------

def generated_antimony(n: int) -> str:
    """The Antimony PyAntiGen generates for an n x n lattice (reactions, rules, parameter values, initial amounts)."""
    import antimony_modules.Histone.methylation_Zheng_PNAS2012 as mod
    from pyantigen.generate.RxnDict_to_antimony import generate_antimony_from_txt
    from pyantigen.generate.pyantigen import PyAntiGen

    mod.N_ME, mod.STATES = n, mod.lattice(n)
    model = PyAntiGen(name=MODEL, isotopes=[""])
    mod.MethylationGrid_Zheng_PNAS2012(model)
    with tempfile.TemporaryDirectory() as tmp:
        f = Path(tmp) / "reaction_dict.txt"
        f.write_text("".join(str(r) + "\n" for r in model.reactions), encoding="utf-8")
        script, species, parameters, _comps, errors = generate_antimony_from_txt(str(f), MODEL)
    assert not errors, errors
    values = [f"{p} = 1" for p in sorted(parameters) if p not in ("inflow",)] + ["inflowp = 0.03", "dilution = 1"]
    values += [f"{s} = {1 / len(species):.6g}" for s in sorted(species)]
    return "\n".join([script.rstrip(), *model.rules, *values]) + "\n"


def exported_sbml(antimony_text: str) -> str:
    import tellurium as te
    return te.antimonyToSBML(antimony_text)


# --------------------------------------------------------------------------------------------------------------------
# the tables
# --------------------------------------------------------------------------------------------------------------------

def fmt(s: Size) -> str:
    return f"{s.bytes:,} | {s.lines:,} | {s.code:,} | {s.gzip:,}"


def x(a: Size, b: Size, field: str) -> str:
    """b over a, as a multiple."""
    return f"{getattr(b, field) / max(getattr(a, field), 1):.1f}x"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--petab", default=str(ROOT.parent / "petab-benchmark-problems" / "problems" / MODEL / "v1"))
    args = ap.parse_args(argv)
    petab = load_petab(args.petab)
    pf = petab.files
    L = []

    # ---------------------------------------------------------------- A. the model
    a_end = class_range(MODULE, "MethylationGrid_Zheng_PNAS2012")[1]
    obs_first, obs_last = class_range(MODULE, "Observables_Zheng_PNAS2012")
    module_model = measure(MODULE, 1, a_end)                               # docstring, lattice helpers, the reaction loops
    module_obs = measure(MODULE, obs_first, obs_last)
    values = measure(AM / f"{MODEL}_parameters.csv") + measure(AM / f"{MODEL}_InitialConditions.csv")

    sbml = measure(petab.sbml_path)
    ant = sbml_to_antimony(petab.sbml_path)
    tell_raw, tell_clean = measure_text(ant.raw), measure_text(ant.text)
    gen4 = generated_antimony(4)
    gen = measure_text(gen4)
    sbml_from_gen = measure_text(exported_sbml(gen4))

    L += [f"# How much does PyAntiGen compress {MODEL}?", "",
          "Two questions, measured separately. Sizes: bytes, non-blank lines, **code lines** (no comments or docstrings; "
          "for data files, the same as lines) and **gzip -9** bytes, the baseline any representation has to beat to claim it "
          "removed information and not only repetition.", "",
          "## A. The model: how much text defines the reactions?", "",
          "| representation | bytes | lines | code | gzip |", "| --- | ---: | ---: | ---: | ---: |",
          f"| PEtab SBML (`model.xml`, one function definition per reaction) | {fmt(sbml)} |",
          f"| SBML exported from PyAntiGen's Antimony (for scale) | {fmt(sbml_from_gen)} |",
          f"| Antimony from Tellurium (`sbmlToAntimony`) | {fmt(tell_raw)} |",
          f"| ... without the `X is \"name\"` lines, annotations and notes | {fmt(tell_clean)} |",
          f"| Antimony PyAntiGen generates (reactions, rules, parameter and initial values) | {fmt(gen)} |",
          f"| **PyAntiGen module: lattice, reaction loops** (`methylation_Zheng_PNAS2012.py`) | {fmt(module_model)} |",
          f"| module + the value files it needs (`_parameters.csv`, `_InitialConditions.csv`) | {fmt(module_model + values)} |",
          "",
          "Compression of the *definition* (rows against the module alone, ratio of that representation to the module):", "",
          "| from | bytes | code lines | gzip |", "| --- | ---: | ---: | ---: |"]
    for name, s in (("PEtab SBML", sbml), ("Tellurium Antimony", tell_raw), ("clean Antimony", tell_clean),
                    ("PyAntiGen's generated Antimony", gen)):
        # lines of pretty-printed XML are formatting, not statements: no line ratio for the SBML
        lines = "n/a" if name == "PEtab SBML" else x(module_model, s, "code")
        L.append(f"| {name} | {x(module_model, s, 'bytes')} | {lines} | {x(module_model, s, 'gzip')} |")
    L += ["", "Reading it honestly: the Antimony PyAntiGen generates is **not** shorter than Tellurium's, it is the same "
          "network written out (`k * (S/V) * V` with the volume explicit). The gain is in what a person has to write and read: "
          "the generator, not its output. gzip shows how much of the SBML is mere repetition: it shrinks "
          f"{sbml.bytes / sbml.gzip:.0f}x on its own, so the module's {sbml.gzip / module_model.gzip:.1f}x advantage over the "
          "gzipped SBML is the part that is more than repetition. Against the generated Antimony the gzip ratio is "
          f"{gen.gzip / module_model.gzip:.1f}x: that file is itself very regular, and the module is little more than its "
          "irreducible content. (Bytes are counted with LF line endings, so they do not depend on the platform; files on a "
          "Windows checkout are a little larger.)", ""]

    # scaling
    import antimony_modules.Histone.methylation_Zheng_PNAS2012 as mod
    L += ["### Growth with the size of the system", "",
          "The same module, `N_ME` changed from 4 to other lattice sizes (an n x n lattice of methylation states): its text "
          "does not change; the Antimony and SBML it stands for grow with the number of reactions.", "",
          "| lattice | states | reactions | module code lines | generated Antimony (bytes) | SBML from it (bytes) | SBML / module |",
          "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    scaling = []
    for n in (2, 3, 4, 6, 8, 12):
        text = generated_antimony(n)
        sb = measure_text(exported_sbml(text))
        g = measure_text(text)
        n_states = n * n - 1
        n_rxn = 2 * (2 * n * (n - 1) - 2) + 1 + n_states               # reversible steps x2, inflow, dilutions
        scaling.append((n, n_states, n_rxn, g.bytes, sb.bytes))
        L.append(f"| {n} x {n} | {n_states} | {n_rxn} | {module_model.code} | {g.bytes:,} | {sb.bytes:,} | "
                 f"{sb.bytes / module_model.bytes:.0f}x |")
    mod.N_ME, mod.STATES = 4, mod.lattice(4)
    L += ["", "(Reactions counted as the engine writes them: each reversible step is two.)", ""]

    # ---------------------------------------------------------------- B. the estimation problem
    pt = {k: measure(pf[k][0]) for k in ("parameters", "conditions", "measurements", "observables", "yaml")}
    viz = measure(pf["visualization"][0])
    py = {n: measure(MODS / f"{n}.py") for n in ("Optimizer_settings", "Loss_config", "Events", "Solver_settings", "Data",
                                                  "Experiment", "Plots", "Update_parameters", "Update_opt_parameters",
                                                  "Observed_species")}
    datacsv = measure(ROOT / "data" / f"{MODEL}.csv")
    params_csv = measure(AM / f"{MODEL}_parameters.csv")
    harness = total([measure(PROJECT / f) for f in ("AntiGen_paths.py", "Model_run.py", "Model_generate.py")])

    pairs = [
        ("parameters (nominal values, bounds, scale, estimate)", pt["parameters"],
         py["Optimizer_settings"] + params_csv, "`Optimizer_settings.py` + `_parameters.csv`"),
        ("observables and noise model", pt["observables"], py["Loss_config"] + module_obs,
         "`Loss_config.py` + the `Observables` class"),
        ("conditions, pre-equilibration", pt["conditions"], py["Events"] + py["Solver_settings"],
         "`Events.py` + `Solver_settings.py`"),
        ("measurements", pt["measurements"], datacsv + py["Data"], "`data/*.csv` + `Data.py`"),
        ("problem definition (what is simulated, what is fitted)", pt["yaml"],
         py["Experiment"] + py["Observed_species"] + py["Update_parameters"] + py["Update_opt_parameters"],
         "`Experiment.py` + three small hooks"),
        ("visualizations", viz, py["Plots"], "`Plots.py`"),
    ]
    petab_total = total(p for _n, p, _q, _f in pairs)
    pa_total = total(q for _n, _p, q, _f in pairs)
    L += ["## B. The estimation problem: PEtab tables against Experiment, Loss_config and the rest", "",
          "Everything PEtab states besides the SBML, against the PyAntiGen files that say the same. The PyAntiGen side also "
          "contains things PEtab has no place for (three optimizer specs, solver and backend settings, a known-sigma loss "
          "variant for Levenberg-Marquardt) and the observables' *formulas*, which PEtab keeps in `observables.tsv`.", "",
          "| PEtab | bytes | lines | gzip | PyAntiGen | bytes | code lines | gzip | bytes ratio |",
          "| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: |"]
    for name, p, q, files in pairs:
        L.append(f"| {name} | {p.bytes:,} | {p.lines:,} | {p.gzip:,} | {files} | {q.bytes:,} | {q.code:,} | {q.gzip:,} | "
                 f"{q.bytes / max(p.bytes, 1):.2f} |")
    L += [f"| **total** | **{petab_total.bytes:,}** | **{petab_total.lines:,}** | **{petab_total.gzip:,}** | | "
          f"**{pa_total.bytes:,}** | **{pa_total.code:,}** | **{pa_total.gzip:,}** | "
          f"**{pa_total.bytes / petab_total.bytes:.2f}** |", "",
          f"On top of that, PyAntiGen needs a fixed harness (`AntiGen_paths.py`, `Model_run.py`, `Model_generate.py`: "
          f"{harness.bytes:,} bytes, {harness.code} code lines) that is the same for every model and has no PEtab counterpart.", "",
          "### What is information and what is boilerplate", "",
          "A good deal of PyAntiGen's problem statement is structure that a model does not choose: the `Optimization` and "
          "`Experiment` dataclasses, the hooks that do nothing here (`Update_parameters`, `Update_opt_parameters`, "
          "`Observed_species`), and the plot code. The lines that *state the problem*:", ""]
    def fn(path, name):
        lo, hi = function_range(path, name)
        return measure(path, lo, hi).code

    info = [
        ("the 15 observables and their noise model", "`Loss_config.py`: `Zheng_loss_config`",
         fn(MODS / "Loss_config.py", "Zheng_loss_config")),
        ("the observables' formulas (fractions of total H3)", "the `Observables` class", module_obs.code),
        ("which parameters are fitted, from where, on what scale, within what bounds",
         "`Optimizer_settings.py`: `_build_Zheng_fit`", fn(MODS / "Optimizer_settings.py", "_build_Zheng_fit")),
        ("the experiment: pre-equilibrate, switch dilution off, observe 25 h",
         "`Solver_settings.py`: `solver_settings_Zheng`; `Events.py`: `Zheng_event`",
         fn(MODS / "Solver_settings.py", "solver_settings_Zheng") + fn(MODS / "Events.py", "Zheng_event")),
        ("which data belong to the experiment", "`Data.py`: `load_Zheng_data`", fn(MODS / "Data.py", "load_Zheng_data")),
    ]
    L += ["| statement | where | code lines |", "| --- | --- | ---: |"] + [f"| {a} | {b} | {c} |" for a, b, c in info]
    state_total = sum(c for _a, _b, c in info)
    L += ["", f"Together {state_total} code lines, against {petab_total.lines:,} lines in PEtab's tables. The rest of "
          f"PyAntiGen's {pa_total.code} code lines is framework (the `Optimization` and `Experiment` dataclasses, the solver "
          "factory, backend settings, plot code, empty hooks) and things PEtab cannot say at all (optimizer choice, "
          "solver tolerances, a second loss for Levenberg-Marquardt).", ""]

    # a closing line that does not oversell
    ratio_b = pa_total.bytes / petab_total.bytes
    L += ["## Summary", "",
          f"* **Model.** The module is {sbml.bytes / module_model.bytes:.0f}x smaller than PEtab's SBML in bytes "
          f"({sbml.gzip / module_model.gzip:.0f}x after gzip, which already removes the SBML's repetition) and "
          f"{gen.bytes / module_model.bytes:.1f}x smaller than the Antimony it generates ({gen.gzip / module_model.gzip:.1f}x "
          "after gzip). The generated Antimony is as long as Tellurium's: nothing is gained in the Antimony itself, only in "
          "the text a person writes. The module does not change as the system grows while the SBML grows with the reactions "
          f"({scaling[-1][4] / module_model.bytes:.0f}x the module at a 12 x 12 lattice).",
          f"* **Estimation problem.** In bytes it is {ratio_b:.2f}x PEtab's tables ({pa_total.bytes:,} against "
          f"{petab_total.bytes:,}): no real saving, because PEtab's tables are compact for one condition and 15 observables "
          f"and PyAntiGen's problem files carry framework structure. The statements themselves are {state_total} code lines "
          f"against {petab_total.lines} table lines: the 15 observables are one loop, the fitted parameters come from the "
          "model's own rate-constant helper, and a different sigma model or optimizer is a one-key change. Conditions are where "
          "PEtab is clearly more compact (one 71-byte table against the solver blocks and an event).", ""]
    (HERE / "COMPRESSION.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    figure(sbml, tell_raw, tell_clean, gen, module_model, scaling, pairs)


# --------------------------------------------------------------------------------------------------------------------
# figure
# --------------------------------------------------------------------------------------------------------------------

def figure(sbml, tell_raw, tell_clean, gen, module_model, scaling, pairs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    PETAB, PYANTIGEN = "#eb6834", "#2a78d6"                 # categorical slots 2 and 1 of the reference palette
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0), facecolor=SURFACE,
                             gridspec_kw={"width_ratios": [1.15, 1.0, 1.25]})
    for ax in axes:
        ax.set_facecolor(SURFACE)
        ax.grid(True, axis="x", color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        ax.tick_params(colors=MUTED, labelsize=8.5)
        for sp in ax.spines.values():
            sp.set_visible(False)

    # A1: the model, bytes
    ax = axes[0]
    rows = [("PEtab SBML", sbml.bytes, PETAB), ("Tellurium Antimony", tell_raw.bytes, PETAB),
            ("clean Antimony", tell_clean.bytes, PETAB), ("PyAntiGen's generated Antimony", gen.bytes, PYANTIGEN),
            ("PyAntiGen module", module_model.bytes, PYANTIGEN)]
    ys = range(len(rows))[::-1]
    ax.barh(list(ys), [r[1] for r in rows], color=[r[2] for r in rows], height=0.6)
    for y, r in zip(ys, rows):
        ax.text(r[1] + sbml.bytes * 0.01, y, f"{r[1]:,} B", va="center", fontsize=8.5, color=INK)
    ax.set_yticks(list(ys))
    ax.set_yticklabels([r[0] for r in rows], color=INK, fontsize=9)
    ax.set_xlim(0, sbml.bytes * 1.22)
    ax.set_title("A. The model: bytes that define the reactions", loc="left", color=INK, fontsize=10.5)

    # A2: growth
    ax = axes[1]
    ns = [s[2] for s in scaling]
    ax.grid(True, axis="y", color=GRID, linewidth=0.8)
    ax.plot(ns, [s[4] for s in scaling], color=PETAB, marker="D", markersize=5, linewidth=2, label="SBML")
    ax.plot(ns, [s[3] for s in scaling], color=PYANTIGEN, marker="o", markersize=5, linewidth=2, alpha=0.55,
            label="PyAntiGen's generated Antimony")
    ax.plot(ns, [module_model.bytes] * len(ns), color=PYANTIGEN, linewidth=2, label="PyAntiGen module")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("reactions in the system", color=MUTED)
    ax.set_ylabel("bytes (log scale)", color=MUTED)
    ax.annotate("SBML", (ns[-1], scaling[-1][4]), textcoords="offset points", xytext=(-6, 6), ha="right", color=INK, fontsize=9)
    ax.annotate("generated Antimony", (ns[-1], scaling[-1][3]), textcoords="offset points", xytext=(-6, 6), ha="right",
                color=INK, fontsize=9)
    ax.annotate("the module (unchanged)", (ns[0], module_model.bytes), textcoords="offset points", xytext=(4, 6),
                color=INK, fontsize=9)
    ax.set_title("The same module as the system grows (2x2 to 12x12)", loc="left", color=INK, fontsize=10.5)

    # B: estimation problem, paired bars per component
    ax = axes[2]
    names = ["parameters", "observables + noise", "conditions", "measurements", "problem definition", "visualizations"]
    k = len(pairs)
    ys = list(range(k))[::-1]
    h = 0.36
    ax.barh([y + h / 2 for y in ys], [p.bytes for _n, p, _q, _f in pairs], height=h, color=PETAB)
    ax.barh([y - h / 2 for y in ys], [q.bytes for _n, _p, q, _f in pairs], height=h, color=PYANTIGEN)
    top = max(max(p.bytes, q.bytes) for _n, p, q, _f in pairs)
    for y, (_n, p, q, _f) in zip(ys, pairs):
        ax.text(p.bytes + top * 0.01, y + h / 2, f"{p.bytes:,}", va="center", fontsize=8, color=INK)
        ax.text(q.bytes + top * 0.01, y - h / 2, f"{q.bytes:,}", va="center", fontsize=8, color=INK)
    ax.set_yticks(ys)
    ax.set_yticklabels(names, color=INK, fontsize=9)
    ax.set_xlim(0, top * 1.15)
    ax.set_title("B. The estimation problem: bytes per component", loc="left", color=INK, fontsize=10.5)
    ax.text(top * 1.13, ys[-1] - 0.55, "orange = PEtab table, blue = PyAntiGen files", ha="right", fontsize=8.5, color=MUTED)

    fig.tight_layout()
    fig.savefig(HERE / "compression.png", dpi=150, facecolor=SURFACE)
    print("wrote", HERE / "compression.png")


if __name__ == "__main__":
    main()

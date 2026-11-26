"""How much does PyAntiGen compress Zhao_QuantBiol2020?  Two questions, measured separately.

  A. WHAT IS SIMULATED.  PEtab: the SBML (3 reactions, plus the rules and initial assignments that pick each region's and
     stage's values with six 0/1 "trigger" parameters) and conditions.tsv.  PyAntiGen: a 3-reaction module, a 7-row
     conditions table, and two small hooks that give each replicate its values.
  B. THE ESTIMATION PROBLEM.  PEtab's parameters / observables / measurements / problem / visualization tables against
     PyAntiGen's Optimizer_settings / Loss_config / Data / Experiment / Plots.

    python benchmarks/Zhao_QuantBiol2020/compression.py [--petab <dir with problem.yaml>]

Writes COMPRESSION.md and compression.png next to this file. Sizes: bytes, non-blank lines, code lines (no comments or
docstrings) and gzip -9 (petab2pyantigen.compression).
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(ROOT))

from petab2pyantigen.compression import Size, function_range, measure, measure_text, total        # noqa: E402
from petab2pyantigen.petab_reader import load_petab                                               # noqa: E402
from petab2pyantigen.sbml_to_antimony import sbml_to_antimony                                     # noqa: E402

MODEL = "Zhao_QuantBiol2020"
MODULE = ROOT / "antimony_modules" / "Epidemic" / "siqc_Zhao_QuantBiol2020.py"
PROJECT = ROOT / "problems" / MODEL
MODS = PROJECT / "Modules"
AM = ROOT / "antimony_models" / MODEL
DATA = ROOT / "data"


def fmt(s: Size) -> str:
    return f"{s.bytes:,} | {s.lines:,} | {s.code:,} | {s.gzip:,}"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--petab", default=str(ROOT.parent / "petab-benchmark-problems" / "problems" / MODEL / "v1"))
    args = ap.parse_args(argv)
    petab = load_petab(args.petab)
    pf = petab.files
    L = []

    # ------------------------------------------------------------------ A. what is simulated
    sbml = measure(petab.sbml_path)
    cond_tsv = measure(pf["conditions"][0])
    ant = sbml_to_antimony(petab.sbml_path)
    tell_raw, tell_clean = measure_text(ant.raw), measure_text(ant.text)

    from pyantigen.generate.antimony_utils import load_antimony_files
    gen = measure_text(load_antimony_files(MODEL, str(ROOT)))               # the Antimony the engine simulates
    module = measure(MODULE)
    values = measure(AM / f"{MODEL}_parameters.csv") + measure(AM / f"{MODEL}_InitialConditions.csv")
    cond_csv = measure(DATA / f"{MODEL}_conditions.csv")
    hooks = measure(MODS / "Update_parameters.py") + measure(MODS / "Update_opt_parameters.py")
    exp_fn = measure(MODS / "Experiment.py", *function_range(MODS / "Experiment.py", "_build_experiment"))
    sim_petab = sbml + cond_tsv
    sim_pa = module + values + cond_csv + hooks + exp_fn

    L += [f"# How much does PyAntiGen compress {MODEL}?", "",
          "Two questions, measured separately. Sizes: bytes, non-blank lines, **code lines** (no comments or docstrings; for "
          "data files, the same as lines) and **gzip -9** bytes, the baseline any representation has to beat to claim it "
          "removed information and not only repetition.", "",
          "## A. What is simulated", "",
          "The model is three reactions. PEtab spreads *which region and stage* a simulation is over several places: the six "
          "`Trigger_*` parameters of `conditions.tsv`, rules `R := Trigger_Wuhan * (Trigger_Stage_I * R_Stage_I_Wuhan + ...)` "
          "(also for gamma_1, gamma_2, the population and alpha, beta), and four initial assignments of the same shape. "
          "PyAntiGen has one replicate per (region, stage) that sets its own values.", "",
          "| representation | bytes | lines | code | gzip |", "| --- | ---: | ---: | ---: | ---: |",
          f"| PEtab SBML (`model.xml`) | {fmt(sbml)} |",
          f"| PEtab `conditions.tsv` | {fmt(cond_tsv)} |",
          f"| **PEtab, what is simulated: SBML + conditions** | **{fmt(sim_petab)}** |",
          f"| Antimony from Tellurium (`sbmlToAntimony`) | {fmt(tell_raw)} |",
          f"| ... without the `X is \"name\"` lines, annotations and notes | {fmt(tell_clean)} |",
          f"| Antimony PyAntiGen generates (3 reactions, rules, parameter and initial values) | {fmt(gen)} |",
          f"| PyAntiGen module: the reactions and rules (`siqc_Zhao_QuantBiol2020.py`) | {fmt(module)} |",
          f"| + value files (`_parameters.csv`, `_InitialConditions.csv`) | {fmt(values)} |",
          f"| + conditions table (`data/*_conditions.csv`, 7 rows) | {fmt(cond_csv)} |",
          f"| + hooks (`Update_parameters.py`, `Update_opt_parameters.py`) | {fmt(hooks)} |",
          f"| + `Experiment._build_experiment` (the loop over conditions) | {fmt(exp_fn)} |",
          f"| **PyAntiGen, what is simulated: all of the above** | **{fmt(sim_pa)}** |", "",
          f"In bytes PyAntiGen's definition is {sim_petab.bytes / sim_pa.bytes:.1f}x smaller than PEtab's "
          f"({sim_pa.bytes:,} against {sim_petab.bytes:,}); after gzip {sim_petab.gzip / sim_pa.gzip:.1f}x "
          f"({sim_pa.gzip:,} against {sim_petab.gzip:,}). Against the Antimony Tellurium writes from the SBML "
          f"({tell_clean.bytes:,} bytes clean), the whole PyAntiGen definition is {tell_clean.bytes / sim_pa.bytes:.1f}x "
          "smaller, and that Antimony still carries all 37 parameters and the trigger rules.", ""]

    # ------------------------------------------------------------------ B. the estimation problem
    pt = {k: measure(pf[k][0]) for k in ("parameters", "measurements", "observables", "yaml")}
    viz = measure(pf["visualization"][0])
    py = {n: measure(MODS / f"{n}.py") for n in ("Optimizer_settings", "Loss_config", "Data", "Experiment", "Plots",
                                                  "Solver_settings", "Events", "Observed_species")}
    datacsv = measure(DATA / f"{MODEL}.csv")
    harness = total([measure(PROJECT / f) for f in ("AntiGen_paths.py", "Model_run.py", "Model_generate.py")])
    pairs = [
        ("parameters (28 estimated: 21 rates, 7 noise; bounds, scale, nominal values)", pt["parameters"],
         py["Optimizer_settings"], "`Optimizer_settings.py` (nominal values are in the conditions table)"),
        ("observables and noise model (5 observables, 1 measured)", pt["observables"], py["Loss_config"],
         "`Loss_config.py` (one observable; one sigma per replicate is the default)"),
        ("measurements", pt["measurements"], datacsv + py["Data"], "`data/*.csv` + `Data.py`"),
        ("problem definition and solver", pt["yaml"], py["Experiment"] + py["Solver_settings"] + py["Events"] + py["Observed_species"],
         "`Experiment.py` + `Solver_settings.py` + two stubs"),
        ("visualizations", viz, py["Plots"], "`Plots.py`"),
    ]
    petab_total = total(p for _n, p, _q, _f in pairs)
    pa_total = total(q for _n, _p, q, _f in pairs)
    L += ["## B. The estimation problem: PEtab tables against Loss_config, Optimizer_settings and the rest", "",
          "Everything PEtab states besides the SBML and the conditions (which part A counts). `Experiment.py` here also holds "
          "the framework's dataclass; the PyAntiGen side contains things PEtab has no place for (the optimizer choice, an "
          "AMICI variant, solver tolerances).", "",
          "| PEtab | bytes | lines | gzip | PyAntiGen | bytes | code lines | gzip | bytes ratio |",
          "| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: |"]
    for name, p, q, files in pairs:
        L.append(f"| {name} | {p.bytes:,} | {p.lines:,} | {p.gzip:,} | {files} | {q.bytes:,} | {q.code:,} | {q.gzip:,} | "
                 f"{q.bytes / max(p.bytes, 1):.2f} |")
    L += [f"| **total** | **{petab_total.bytes:,}** | **{petab_total.lines:,}** | **{petab_total.gzip:,}** | | "
          f"**{pa_total.bytes:,}** | **{pa_total.code:,}** | **{pa_total.gzip:,}** | **{pa_total.bytes / petab_total.bytes:.2f}** |", "",
          f"A fixed harness (`AntiGen_paths.py`, `Model_run.py`, `Model_generate.py`: {harness.bytes:,} bytes, {harness.code} code "
          "lines), the same for every model, has no PEtab counterpart.", ""]

    def fn(path, name):
        lo, hi = function_range(path, name)
        return measure(path, lo, hi).code

    info = [
        ("the observable and its noise", "`Loss_config.py`: `Zhao_loss_config`", fn(MODS / "Loss_config.py", "Zhao_loss_config")),
        ("the 21 fitted parameters, their bounds, start and scale", "`Optimizer_settings.py`: `_build_Zhao_fit`",
         fn(MODS / "Optimizer_settings.py", "_build_Zhao_fit")),
        ("what distinguishes the 7 conditions, and how a fitted name reaches the model",
         "`Update_parameters.py`: `update_Zhao`; `Update_opt_parameters.py`: `update_opt_Zhao`",
         fn(MODS / "Update_parameters.py", "update_Zhao") + fn(MODS / "Update_opt_parameters.py", "update_opt_Zhao")),
        ("one replicate per condition", "`Experiment.py`: `_build_experiment`", fn(MODS / "Experiment.py", "_build_experiment")),
        ("which data belong to a replicate", "`Data.py`: `load_Zhao_data`", fn(MODS / "Data.py", "load_Zhao_data")),
        ("the simulated days", "`Solver_settings.py`: `solver_settings_Zhao`", fn(MODS / "Solver_settings.py", "solver_settings_Zhao")),
    ]
    state_total = sum(c for _a, _b, c in info)
    L += ["### The statements themselves", "", "| statement | where | code lines |", "| --- | --- | ---: |"]
    L += [f"| {a} | {b} | {c} |" for a, b, c in info]
    L += ["", f"Together {state_total} code lines, plus the reaction module ({module.code} code lines) and two small data files "
          f"({cond_csv.lines} + {datacsv.lines} lines).", ""]

    all_petab = petab_total + sim_petab
    all_pa = pa_total + sim_pa
    L += ["## Summary", "",
          f"* **What is simulated:** {sim_petab.bytes / sim_pa.bytes:.1f}x smaller in bytes ({sim_petab.gzip / sim_pa.gzip:.1f}x "
          f"after gzip). PEtab's SBML is {sbml.bytes / module.bytes:.0f}x the reaction module because it also carries the "
          "trigger rules and initial assignments that PyAntiGen replaces with a seven-row table and two hooks.",
          f"* **Estimation problem:** {pa_total.bytes / petab_total.bytes:.2f}x the bytes of PEtab's tables "
          f"({pa_total.bytes:,} against {petab_total.bytes:,}), and almost all of the saving is the measurement table: PEtab "
          f"repeats the observable, the condition and the noise parameter on each of its 82 rows ({pt['measurements'].bytes:,} "
          f"bytes against {datacsv.bytes + py['Data'].bytes:,} for the data file and its loader). The parameters and the "
          "problem definition are larger in PyAntiGen, which carries framework structure there.",
          f"* **Everything, harness excluded:** {all_pa.bytes:,} bytes against {all_petab.bytes:,} "
          f"({all_petab.bytes / all_pa.bytes:.1f}x).", ""]
    (HERE / "COMPRESSION.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    figure(sbml, cond_tsv, sim_petab, tell_clean, gen, module, sim_pa, pairs)


def figure(sbml, cond_tsv, sim_petab, tell_clean, gen, module, sim_pa, pairs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    PETAB, PYANTIGEN = "#eb6834", "#2a78d6"
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.6), facecolor=SURFACE, gridspec_kw={"width_ratios": [1.1, 1.0]})
    for ax in axes:
        ax.set_facecolor(SURFACE)
        ax.grid(True, axis="x", color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        ax.tick_params(colors=MUTED, labelsize=8.5)
        for sp in ax.spines.values():
            sp.set_visible(False)

    ax = axes[0]
    rows = [("PEtab: SBML + conditions", sim_petab.bytes, PETAB), ("Tellurium Antimony (clean)", tell_clean.bytes, PETAB),
            ("PyAntiGen: generated Antimony", gen.bytes, PYANTIGEN), ("PyAntiGen: everything that defines\nwhat is simulated", sim_pa.bytes, PYANTIGEN),
            ("PyAntiGen: reaction module", module.bytes, PYANTIGEN)]
    ys = list(range(len(rows)))[::-1]
    ax.barh(ys, [r[1] for r in rows], color=[r[2] for r in rows], height=0.6)
    for y, r in zip(ys, rows):
        ax.text(r[1] + sim_petab.bytes * 0.01, y, f"{r[1]:,} B", va="center", fontsize=8.5, color=INK)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], color=INK, fontsize=9)
    ax.set_xlim(0, sim_petab.bytes * 1.2)
    ax.set_title("A. What is simulated: bytes", loc="left", color=INK, fontsize=10.5)

    ax = axes[1]
    names = ["parameters", "observables + noise", "measurements", "problem + solver", "visualizations"]
    ys = list(range(len(pairs)))[::-1]
    h = 0.36
    top = max(max(p.bytes, q.bytes) for _n, p, q, _f in pairs)
    ax.barh([y + h / 2 for y in ys], [p.bytes for _n, p, _q, _f in pairs], height=h, color=PETAB)
    ax.barh([y - h / 2 for y in ys], [q.bytes for _n, _p, q, _f in pairs], height=h, color=PYANTIGEN)
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

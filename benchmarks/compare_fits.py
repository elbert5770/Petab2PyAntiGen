"""Compare the routes that fit Zheng_PNAS2012: PyAntiGen on RoadRunner / AMICI, and PEtab straight through AMICI.

    python benchmarks/compare_fits.py [--petab <dir with problem.yaml>]

Reads benchmarks/fits/*.json (written by fit_zheng.py and amici_direct_zheng.py) and scores every fitted parameter vector
with ONE objective, PEtab's own: the problem simulated from the SBML and tables by RoadRunner with steady-state
pre-equilibration (petab2pyantigen.reference), sigma at its maximum-likelihood value.  So the routes are compared on the
same footing whatever engine produced them, and a route that reports a flattering number of its own cannot hide.

Writes benchmarks/COMPARISON.md and benchmarks/comparison.png.
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from petab2pyantigen.petab_reader import load_petab                          # noqa: E402
from petab2pyantigen.reference import nominal_values, simulate_petab           # noqa: E402

N = 60                                                       # data points in Zheng_PNAS2012
CONST = 0.5 * N * (1 + math.log(2 * math.pi))                # -log L = NLL (engine convention) + CONST at the ML sigma


def score(problem, params: dict) -> dict:
    """PEtab reference objective at ``params``: SSE, ML sigma, NLL in the engine's convention and as -log L."""
    sims = simulate_petab(problem, {k: v for k, v in params.items() if k != "sigma"})
    sse = float(((sims["measurement"] - sims["simulation"]) ** 2).sum())
    nll = 0.5 * N * math.log(sse / N)
    return {"sse": sse, "sigma_ml": math.sqrt(sse / N), "nll": nll, "neg_log_l": nll + CONST}


def parameter_agreement(fit: dict, nominal: dict) -> dict:
    """How far, in decades, the fitted rate constants are from PEtab's published values."""
    names = [k for k in fit["params"] if k != "sigma" and k in nominal]
    d = np.array([abs(math.log10(fit["params"][k]) - math.log10(nominal[k])) for k in names])
    return {"median": float(np.median(d)), "max": float(d.max()), "n_over_0p1": int((d > 0.1).sum()), "n": len(names)}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--petab", default=str(ROOT.parent / "petab-benchmark-problems" / "problems" / "Zheng_PNAS2012" / "v1"))
    args = ap.parse_args(argv)

    problem = load_petab(args.petab)
    nominal = nominal_values(problem)
    base = score(problem, nominal)
    fits = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((HERE / "fits").glob("*.json"))]
    if not fits:
        sys.exit("no fits in benchmarks/fits: run fit_zheng.py and amici_direct_zheng.py first")

    rows = []
    for f in fits:
        sc = score(problem, f["params"])
        ag = parameter_agreement(f, nominal)
        rows.append({**f, **sc, "d_nll": sc["nll"] - base["nll"], **{f"ag_{k}": v for k, v in ag.items()},
                     "s_per_eval": f["wall_s"] / max(f["nfev"] or 1, 1)})
    rows.sort(key=lambda r: r["wall_s"])

    lines = [
        "# Zheng_PNAS2012: the fitting routes compared",
        "",
        "Every route starts from the same perturbed point (seed 1: each parameter x 10**U(-0.3, 0.3)) and fits all 45 rate "
        "constants and `inflowp` on a log10 scale in [1e-5, 1e3]. Each fitted parameter vector is then **re-scored with one "
        "objective, PEtab's own** (SBML + tables, RoadRunner, steady-state pre-equilibration, sigma at its ML value), so the "
        "routes are compared on the same footing.",
        "",
        f"PEtab's published parameters score NLL = {base['nll']:.4f} (-log L = {base['neg_log_l']:.4f}, "
        f"sigma = {base['sigma_ml']:.6f}). NLL is in the engine's convention, `(n/2) ln(SSE/n)`; -log L adds n/2 (1 + ln 2 pi).",
        "",
        "| route | sigma | converged | evals | wall s | s/eval | NLL (PEtab reference) | vs PEtab's | median / max shift in decades | > 0.1 decade |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |",
    ]
    for r in rows:
        lines.append(
            f"| {r['route']}{'' if r.get('exact_derivatives', True) else ' (finite-difference derivatives)'} | {r['sigma']} | "
            f"{'yes' if r['success'] else 'no (budget)'} | {r['nfev']} | {r['wall_s']:.0f} | {r['s_per_eval']:.3f} | "
            f"{r['nll']:.4f} | {r['d_nll']:+.4f} | {r['ag_median']:.3f} / {r['ag_max']:.2f} | {r['ag_n_over_0p1']} of {r['ag_n']} |")
    lines += ["", "Platform / versions: " + "; ".join(sorted({f"{r['platform']}, Python {r['python']}" for r in rows})), ""]
    notes = [r for r in rows if r.get("note")]
    if notes:
        lines += ["Notes: " + " ".join(f"**{r['label']}**: {r['note']}" for r in notes), ""]
    (HERE / "COMPARISON.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    figure(rows, base)


def figure(rows, base):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    PYANTIGEN, DIRECT = "#2a78d6", "#eb6834"                  # categorical slots 1 and 2 of the reference palette
    markers = {"roadrunner": "o", "amici": "s"}
    fig, ax = plt.subplots(figsize=(9.0, 5.0), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    # label placement per route (offset in points, alignment): the converged routes sit close together
    place = {"amici_lm": ((10, -14), "left", "top"), "direct_amici": ((-9, 12), "right", "bottom"),
             "amici_lbfgsb": ((9, 12), "left", "bottom"), "neldermead": ((-10, -6), "right", "top")}
    for r in rows:
        direct = r["route"].startswith("PEtab")
        m = "D" if direct else ("^" if r["method"] == "LM" else markers.get(r["backend"], "o"))
        ax.scatter(r["wall_s"], r["d_nll"], s=90, marker=m, color=DIRECT if direct else PYANTIGEN,
                   edgecolor=SURFACE, linewidth=1.5, zorder=3)
        short = r["route"].replace("PyAntiGen / ", "PyAntiGen, ").replace(" / ", ", ")
        off, ha, va = next((v for k, v in place.items() if k in r["label"]), ((9, 7), "left", "bottom"))
        ax.annotate(f"{short}\n{r['nfev']} evals, {r['wall_s']:.0f} s", (r["wall_s"], r["d_nll"]), textcoords="offset points",
                    xytext=off, fontsize=8.5, color=INK, ha=ha, va=va)
    ax.axhline(0, color=MUTED, linewidth=1)
    ax.set_xlim(0, max(r["wall_s"] for r in rows) * 1.25)
    ax.text(ax.get_xlim()[1], 0, "PEtab's published parameters  ", color=MUTED, fontsize=8, va="bottom", ha="right")
    ax.set_yscale("symlog", linthresh=0.01)
    ax.set_xlabel("wall-clock time of the fit (s)", color=MUTED)
    ax.set_ylabel("NLL above PEtab's published fit\n(PEtab reference objective)", color=MUTED)
    ax.set_title("Zheng_PNAS2012 from one perturbed start: time against the quality of the fit", loc="left", color=INK,
                 fontsize=11)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.tick_params(colors=MUTED)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.margins(y=0.3)
    fig.text(0.01, 0.01, "blue = PyAntiGen (circle RoadRunner, square AMICI, triangle AMICI Levenberg-Marquardt); "
             "orange diamond = PEtab straight through AMICI", fontsize=8, color=MUTED)
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.savefig(HERE / "comparison.png", dpi=150, facecolor=SURFACE)
    print("wrote", HERE / "comparison.png")


if __name__ == "__main__":
    main()

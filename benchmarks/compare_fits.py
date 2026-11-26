"""Compare the routes that fit a converted problem: PyAntiGen on RoadRunner / AMICI, and PEtab straight through AMICI.

    python benchmarks/compare_fits.py --problem Zhao_QuantBiol2020 [--petab <dir with problem.yaml>]

Reads benchmarks/<problem>/fits/*.json (written by fit_problem.py and amici_direct.py) and scores every fitted parameter
vector with ONE objective, PEtab's own: the problem simulated from the SBML and tables by RoadRunner (with its
pre-equilibration, if it has one: petab2pyantigen.reference), each noise parameter at its maximum-likelihood value.  So the
routes are compared on the same footing whatever engine produced them, and a route that reports a flattering number of its
own cannot hide.

Writes benchmarks/<problem>/COMPARISON.md and comparison.png.
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

def noise_names(problem) -> set:
    """The noise parameters of the problem: what its measurement table names in ``noiseParameters``."""
    col = problem.measurements.get("noiseParameters")
    return {str(v) for v in col.dropna().unique()} | {"sigma"} if col is not None else {"sigma"}


def score(problem, params: dict) -> dict:
    """PEtab reference objective at ``params``: SSE, ML sigma, NLL in the engine's convention and as -log L.

    Measurements that share a noise parameter share a sigma, taken at its maximum-likelihood value (Zheng: one for all 60
    points; Zhao: one per dataset). NLL is ``sum_g (n_g/2) ln(SSE_g/n_g)``, the engine's concentrated NLL; -log L adds
    ``N/2 (1 + ln 2 pi)``.
    """
    noise = noise_names(problem)
    sims = simulate_petab(problem, {k: v for k, v in params.items() if k not in noise})
    key = sims["noiseParameters"].astype(str) if "noiseParameters" in sims else "sigma"
    nll, sse_total, sigmas = 0.0, 0.0, []
    for _g, rows in sims.groupby(key):
        sse = float(((rows["measurement"] - rows["simulation"]) ** 2).sum())
        n = len(rows)
        nll += 0.5 * n * math.log(sse / n)
        sse_total += sse
        sigmas.append(math.sqrt(sse / n))
    n_all = len(sims)
    return {"sse": sse_total, "sigma_ml": sigmas[0] if len(sigmas) == 1 else float(np.median(sigmas)), "n_sigmas": len(sigmas),
            "nll": nll, "neg_log_l": nll + 0.5 * n_all * (1 + math.log(2 * math.pi))}


def deriv_note(r: dict) -> str:
    """How the route got its derivatives, when that is not the plain case (records before `derivatives` existed too)."""
    d = r.get("derivatives") or ("chain rule, exact" if r.get("exact_derivatives", True) else "finite differences of the whole loss")
    callable_obs = r.get("observables") == "callable"
    if r["backend"] == "roadrunner" or r["route"].startswith("PEtab"):
        return ""
    if r["method"] == "LM" or callable_obs or "finite" in d:
        return f" ({'callable observables, ' if callable_obs else ''}{d})"
    return ""


def parameter_agreement(fit: dict, nominal: dict, noise: set) -> dict:
    """How far, in decades, the fitted parameters are from PEtab's published values (noise parameters excluded)."""
    names = [k for k in fit["params"] if k not in noise and k in nominal]
    d = np.array([abs(math.log10(fit["params"][k]) - math.log10(nominal[k])) for k in names])
    return {"median": float(np.median(d)), "max": float(d.max()), "n_over_0p1": int((d > 0.1).sum()), "n": len(names)}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--problem", default="Zheng_PNAS2012")
    ap.add_argument("--petab", default=None, help="PEtab directory (default: the sibling petab-benchmark-problems checkout)")
    args = ap.parse_args(argv)
    out = HERE / args.problem
    petab = args.petab or str(ROOT.parent / "petab-benchmark-problems" / "problems" / args.problem / "v1")

    problem = load_petab(petab)
    nominal = nominal_values(problem)
    noise = noise_names(problem)
    base = score(problem, nominal)
    fits = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((out / "fits").glob("*.json"))]
    if not fits:
        sys.exit(f"no fits in {out / 'fits'}: run fit_problem.py and amici_direct.py first")

    rows = []
    for f in fits:
        sc = score(problem, f["params"])
        ag = parameter_agreement(f, nominal, noise)
        rows.append({**f, **sc, "d_nll": sc["nll"] - base["nll"], **{f"ag_{k}": v for k, v in ag.items()},
                     "s_per_eval": f["wall_s"] / max(f["nfev"] or 1, 1)})
    rows.sort(key=lambda r: r["wall_s"])

    lines = [
        f"# {args.problem}: the fitting routes compared",
        "",
        "Every route starts from the same perturbed point (seed 1: each parameter x 10**U(-0.3, 0.3), inside its bounds) and "
        "fits the problem's parameters on a log10 scale. Each fitted parameter vector is then **re-scored with one objective, "
        "PEtab's own** (SBML + tables, RoadRunner, each noise parameter at its ML value), so the routes are compared on the "
        "same footing.",
        "",
        f"PEtab's published parameters score NLL = {base['nll']:.4f} (-log L = {base['neg_log_l']:.4f}; "
        f"{base['n_sigmas']} noise parameter{'s' if base['n_sigmas'] > 1 else ''} at ML, "
        f"{'sigma' if base['n_sigmas'] == 1 else 'median sigma'} = {base['sigma_ml']:.6g}). NLL is in the engine's convention, "
        "the sum over noise groups of `(n/2) ln(SSE/n)`; -log L adds N/2 (1 + ln 2 pi).",
        "",
        "| route | sigma | converged | evals | wall s | s/eval | NLL (PEtab reference) | vs PEtab's | median / max shift in decades | > 0.1 decade |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |",
    ]
    for r in rows:
        lines.append(
            f"| {r['route']}{deriv_note(r)} | {r['sigma']} | "
            f"{'yes' if r['success'] else 'no (budget)'} | {r['nfev']} | {r['wall_s']:.0f} | {r['s_per_eval']:.3f} | "
            f"{r['nll']:.4f} | {r['d_nll']:+.4f} | {r['ag_median']:.3f} / {r['ag_max']:.2f} | {r['ag_n_over_0p1']} of {r['ag_n']} |")
    lines += ["", "Platform / versions: " + "; ".join(sorted({f"{r['platform']}, Python {r['python']}" for r in rows})), ""]
    notes = [r for r in rows if r.get("note")]
    if notes:
        lines += ["Notes: " + " ".join(f"**{r['label']}**: {r['note']}" for r in notes), ""]
    (out / "COMPARISON.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    figure(rows, base, out, args.problem)


def figure(rows, base, out, problem_name):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    PYANTIGEN, DIRECT = "#2a78d6", "#eb6834"                  # categorical slots 1 and 2 of the reference palette
    rows = sorted(rows, key=lambda r: r["wall_s"])
    fig, ax = plt.subplots(figsize=(11.5, 0.62 * len(rows) + 2.3), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ys = list(range(len(rows)))[::-1]
    top = max(r["wall_s"] for r in rows)
    labels = []
    for y, r in zip(ys, rows):
        direct = r["route"].startswith("PEtab")
        ax.barh(y, r["wall_s"], height=0.62, color=DIRECT if direct else PYANTIGEN,
                hatch=None if r["success"] else "///", edgecolor=SURFACE if r["success"] else "white", linewidth=0)
        gap = f"{r['d_nll']:+.4f}"
        ax.text(r["wall_s"] + top * 0.012, y, f"{r['wall_s']:.0f} s  |  {r['nfev']:,} evals  |  NLL {gap}"
                + ("  (stopped at its evaluation limit)" if not r["success"] else ""),
                va="center", fontsize=9, color=INK)
        short = r["route"].replace("PyAntiGen / ", "PyAntiGen, ").replace("PEtab -> AMICI directly", "PEtab straight to AMICI")
        short = short.replace(" / ", ", ")
        note = deriv_note(r).strip(" ()")
        labels.append(short + (f"\n{note}" if note else ""))
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=9, color=INK)
    ax.set_xlim(0, top * 1.62)
    ax.set_xlabel("wall-clock time of the fit (s)", color=MUTED)
    ax.grid(True, axis="x", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTED, length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(f"{problem_name} from one perturbed start: time to fit, by route", loc="left", color=INK, fontsize=11.5)
    fig.text(0.01, 0.012, "blue = PyAntiGen, orange = PEtab straight through AMICI with no PyAntiGen. NLL is the PEtab reference "
             "objective of the fitted parameters relative to PEtab's published fit:\n0 = as good, negative = better, positive = worse. "
             "Hatched: stopped at its evaluation limit. One start per route: where the landscape has several optima, "
             "routes can end in different ones.", fontsize=8, color=MUTED, va="bottom")
    fig.tight_layout(rect=[0, 0.07, 1, 1])
    fig.savefig(out / "comparison.png", dpi=150, facecolor=SURFACE)
    print("wrote", out / "comparison.png")


if __name__ == "__main__":
    main()

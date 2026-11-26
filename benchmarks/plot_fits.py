"""One figure per problem: every route's FINAL fit against the data, one panel per condition (or observable).

    python benchmarks/plot_fits.py --problem Zhao_QuantBiol2020 [--petab <dir with problem.yaml>]

Reads benchmarks/<problem>/fits/*.json (written by fit_problem.py and amici_direct.py) and re-simulates every fitted
parameter vector with ONE simulator, PEtab's own reference (SBML + tables, RoadRunner), on a dense time grid, so the curves
are comparable whatever engine produced them. PEtab's published parameters are drawn too (grey, dashed). The legend carries
each route's NLL from the same reference, so a curve that looks poor can be matched to its number.

Writes benchmarks/<problem>/fits.png.
"""
import argparse
import dataclasses
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))

from petab2pyantigen.petab_reader import load_petab                          # noqa: E402
from petab2pyantigen.reference import nominal_values, simulate_petab           # noqa: E402
from compare_fits import deriv_note, noise_names, score                                    # noqa: E402

INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
# A route keeps its colour in every problem (colour follows the entity, not its rank): blue / orange as in comparison.png.
# (a route keeps its colour at any start: the seed-2 labels share the colour of the route)
ROUTE_COLOR = {"pyantigen_amici_lbfgsb": "#2a78d6", "pyantigen_amici_lbfgsb_seed2": "#2a78d6",
               "direct_amici_lbfgsb": "#eb6834", "direct_amici_lbfgsb_seed2": "#eb6834",
               "pyantigen_roadrunner_neldermead": "#1baf7a",
               "pyantigen_amici_lm": "#6250d6", "pyantigen_amici_lm_seed2": "#6250d6",
               "pyantigen_amici_lm_finite_difference": "#e87ba4",
               "pyantigen_amici_lm_callable_observable": "#eda100", "pyantigen_amici_lm_callable_finite_difference": "#e34948"}
SPARE = ["#008300", "#52514e", "#8a5a2b"]


def dense(problem, params: dict, per_unit: float):
    """Reference simulation on a grid finer than the measurements: one row per (condition, observable, time)."""
    m = problem.measurements
    extra = []
    for (pre, sim, obs), g in m.groupby(["preequilibrationConditionId", "simulationConditionId", "observableId"],
                                        dropna=False, sort=False):
        lo, hi = g["time"].min(), g["time"].max()
        n = int(round((hi - lo) * per_unit)) + 1
        for t in np.linspace(lo, hi, n):
            extra.append({**g.iloc[0].to_dict(), "time": float(t), "measurement": np.nan})
    fine = pd.DataFrame(extra)
    noise = noise_names(problem)
    sims = simulate_petab(dataclasses.replace(problem, measurements=fine), {k: v for k, v in params.items() if k not in noise})
    return sims


def panel_keys(problem):
    """Panels: one per simulation condition when each condition has one observable (Zhao), one per observable when there is one
    condition (Zheng), and one per (condition, observable) otherwise (Weber: observables of very different scales)."""
    m = problem.measurements
    conds = list(dict.fromkeys(m["simulationConditionId"]))
    pairs = list(dict.fromkeys(zip(m["simulationConditionId"], m["observableId"])))
    if len(conds) > 1 and len(pairs) == len(conds):
        return [(c, None) for c in conds]
    return pairs


def short(route: str) -> str:
    s = route.replace("PyAntiGen / ", "PyAntiGen, ").replace("PEtab -> AMICI directly", "PEtab straight to AMICI")
    return s.replace(" / ", ", ")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--problem", default="Zheng_PNAS2012")
    ap.add_argument("--petab", default=None)
    args = ap.parse_args(argv)
    out = HERE / args.problem
    petab = args.petab or str(ROOT.parent / "petab-benchmark-problems" / "problems" / args.problem / "v1")

    problem = load_petab(petab)
    fits = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((out / "fits").glob("*.json"))]
    if not fits:
        sys.exit(f"no fits in {out / 'fits'}")

    span = problem.measurements["time"].max() - problem.measurements["time"].min()
    per_unit = max(1.0, 100.0 / span)                     # about 100 points per curve, at least one per time unit
    curves = {"published": dense(problem, nominal_values(problem), per_unit)}
    nll = {"published": score(problem, nominal_values(problem))["nll"]}
    for f in fits:
        try:
            curves[f["label"]] = dense(problem, f["params"], per_unit)
            nll[f["label"]] = score(problem, f["params"])["nll"]
        except RuntimeError:                      # the reference solver cannot integrate this fitted vector (stiff corner of the box)
            curves[f["label"]], nll[f["label"]] = None, math.nan
            print(f"  {f['label']}: the reference solver cannot simulate this parameter vector; not drawn")

    labels = {"published": "PEtab's published parameters"}
    several_starts = len({f.get("seed", 1) for f in fits}) > 1
    labels.update({f["label"]: short(f["route"]) + deriv_note(f) + (f" [start {f.get('seed', 1)}]" if several_starts else "")
                   for f in fits})
    drawn = [k for k in curves if k != "published" and curves[k] is not None]
    failed = [k for k in curves if k != "published" and curves[k] is None]
    order = sorted(drawn, key=lambda k: nll[k])
    best = min(nll[k] for k in order)
    spare = iter(SPARE)
    colors = {k: ROUTE_COLOR.get(k) or next(spare) for k in order}

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panels = panel_keys(problem)
    ncols = 4
    spread = max(nll[k] for k in order) - best
    dec = 3 if spread < 5 else 1                           # show the digits that tell the routes apart
    long_legend = len(order) > 4                           # many routes: the legend goes under the panels, in two columns
    nrows = math.ceil((len(panels) + (0 if long_legend else 1)) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, 2.7 * nrows + 0.8), facecolor=SURFACE, squeeze=False)
    axes = axes.ravel()
    meas = problem.measurements
    conds_all = list(dict.fromkeys(meas["simulationConditionId"]))
    for ax, (cond, obs) in zip(axes, panels):
        ax.set_facecolor(SURFACE)
        sel = (meas["simulationConditionId"] == cond) & ((meas["observableId"] == obs) if obs else True)
        d = meas[sel]
        for key in ["published"] + order[::-1]:           # worst route first, widest: the best is the thin line on top
            c = curves[key]
            cs = c[(c["simulationConditionId"] == cond) & ((c["observableId"] == obs) if obs else True)]
            for _o, g in cs.groupby("observableId"):
                g = g.sort_values("time")
                mark = dict(marker="o", markersize=7, markerfacecolor="none") if len(g) == 1 else {}      # one time point: no line to see
                if key == "published":
                    ax.plot(g["time"], g["simulation"], color="#8a8985", lw=1.3, ls=(0, (4, 3)), zorder=2, **mark)
                else:
                    rank = order.index(key)
                    ax.plot(g["time"], g["simulation"], color=colors[key], zorder=3,
                            lw=1.4 + 1.8 * rank / max(len(order) - 1, 1), **mark)
        ax.scatter(d["time"], d["measurement"], s=16, color=INK, zorder=5, edgecolor=SURFACE, linewidth=0.6)
        title = ((obs.replace("observable_", "") if len(conds_all) == 1 else f"{cond.replace('model1_', '')}: {obs}") if obs
                 else cond.replace("model1_", "").replace("_", " "))
        ax.set_title(title, loc="left", fontsize=9.5, color=INK)
        ax.grid(True, color=GRID, linewidth=0.7)
        ax.set_axisbelow(True)
        ax.tick_params(colors=MUTED, labelsize=8, length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
    for ax in axes[len(panels):]:
        ax.axis("off")
    # legend goes in the first free cell
    leg = axes[len(panels)] if not long_legend and len(panels) < len(axes) else None
    handles = [plt.Line2D([], [], color="#8a8985", lw=1.3, ls=(0, (4, 3)),
                          label=f"{labels['published']}  (NLL {nll['published']:+.{dec}f})"),
               plt.Line2D([], [], marker="o", color=INK, lw=0, markersize=5, label="measurements")]
    handles += [plt.Line2D([], [], color=colors[k], lw=2.2, label=f"{labels[k]}  (NLL {nll[k]:+.{dec}f})") for k in order]
    for k in failed:                                  # drawn as a legend entry only: the reference solver cannot simulate it
        handles.append(plt.Line2D([], [], color="none", label=f"{labels[k]}: not drawn, the reference solver fails there"))
    if leg is not None:
        leg.legend(handles=handles, loc="center left", frameon=False, fontsize=7.6, labelcolor=INK, handlelength=2.2)
    else:
        fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.03), frameon=False, fontsize=8, ncol=2,
                   labelcolor=INK)
    for ax in axes[:len(panels)]:
        ax.set_xlabel("time", color=MUTED, fontsize=8, labelpad=2)
    fig.suptitle(f"{args.problem}: final fit of every route, re-simulated with PEtab's reference simulator",
                 x=0.01, ha="left", fontsize=12, color=INK)
    fig.text(0.01, 0.005, f"Lower NLL is better; routes within {spread:.2f} NLL of each other overlay. One start per curve, "
             "so routes may end in different optima." if spread < 5 else
             "Lower NLL is better. One start per curve: routes end in different optima, and the published parameters were "
             "not fitted with this objective's free noise parameters.", fontsize=8, color=MUTED, va="bottom")
    fig.tight_layout(rect=[0, 0.035 + (0.11 if long_legend else 0), 1, 0.97])
    fig.savefig(out / "fits.png", dpi=150, facecolor=SURFACE)
    print("wrote", out / "fits.png")
    for k in order:
        print(f"  {nll[k]:+10.4f}  {labels[k]}")
    print(f"  {nll['published']:+10.4f}  {labels['published']}")


if __name__ == "__main__":
    main()

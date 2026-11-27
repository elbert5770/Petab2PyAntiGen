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
import re
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
# the smoothed-switch variant of a route (Smith): the same hue, darker
ROUTE_COLOR.update({"pyantigen_amici_lm_smooth": "#3b2a9e", "pyantigen_amici_lm_smooth_seed2": "#3b2a9e",
                    "pyantigen_amici_lbfgsb_smooth": "#14478a", "pyantigen_amici_lbfgsb_smooth_seed2": "#14478a"})
SPARE = ["#008300", "#52514e", "#8a5a2b"]


def dense(problem, params: dict, per_unit: float):
    """Reference simulation on a grid finer than the measurements: one row per (condition, observable, time)."""
    m = problem.measurements.copy()
    if "preequilibrationConditionId" not in m:               # a problem without pre-equilibration (Smith)
        m["preequilibrationConditionId"] = np.nan
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
    """Panels, each a dict ``conds`` (the conditions drawn), ``obs`` (one observable, or None for all of them), ``xmap`` (None: time on
    the x axis; else the x value of each condition) and ``title``.

    One per simulation condition when each condition has one observable (Zhao), one per observable when there is one condition
    (Zheng), one per (condition, observable) otherwise (Weber: observables of very different scales). With many pairs (Smith: 36
    conditions, mostly measured once) one per observable: a time course if the observable has several times in a condition, else
    its value against the condition quantity (a dose) that varies the most among the conditions that measure it."""
    m = problem.measurements
    conds = list(dict.fromkeys(m["simulationConditionId"]))
    pairs = list(dict.fromkeys(zip(m["simulationConditionId"], m["observableId"])))
    if len(conds) > 1 and len(pairs) == len(conds):
        return [{"conds": [c], "obs": None, "xmap": None, "title": c.replace("model1_", "").replace("_", " ")} for c in conds]
    if len(pairs) <= 20:
        return [{"conds": [c], "obs": o, "xmap": None,
                 "title": o.replace("observable_", "") if len(conds) == 1 else f"{c.replace('model1_', '')}: {o}"} for c, o in pairs]
    out = []
    # one observable measured in several series that differ in their observable parameters (Schwen's four ELISA series, each with its own
    # offset, scale and saturation constant): one panel per series
    keys = ["observableId"]
    if "observableParameters" in m and (m.groupby("observableId")["observableParameters"].nunique() > 1).any():
        keys.append("observableParameters")
    for key, g in m.groupby(keys, sort=False):
        obs, series = (key[0], key[1]) if len(keys) == 2 else (key, None)
        used = list(dict.fromkeys(g["simulationConditionId"]))
        suffix = ""
        if series is not None and m[m["observableId"] == obs]["observableParameters"].nunique() > 1:
            suffix = " (" + (re.search(r"nExpID\d", str(series)).group(0) if re.search(r"nExpID\d", str(series)) else str(series)) + ")"
        if g.groupby("simulationConditionId")["time"].nunique().max() > 1:           # a time course (one condition)
            out.append({"conds": used, "obs": obs, "series": series if suffix else None, "xmap": None,
                        "title": obs.replace("__", " ").replace("observable_", "") + suffix})
            continue
        table = problem.conditions.loc[used].apply(pd.to_numeric, errors="coerce")
        column = max(table.columns, key=lambda c: table[c].nunique())
        if len(used) <= 4:                    # a few conditions that differ in several quantities (Smith 3B): one tick each, labelled
            varying = [c for c in table.columns if table[c].nunique() > 1]
            labels = ["\n".join(f"{c.replace('extracellular_ROS', 'H2O2')} {table.loc[u, c]:g}" for c in varying) or u for u in used]
            out.append({"conds": used, "obs": obs, "xmap": {c: i for i, c in enumerate(used)}, "xlabel": "condition",
                        "ticks": labels, "title": obs.replace("__", " ")})
        elif table[column].nunique() > 1:
            out.append({"conds": used, "obs": obs, "xmap": table[column].to_dict(), "xlabel": column,
                        "title": obs.replace("__", " ")})
        else:                                                                       # nothing varies: one tick per condition
            out.append({"conds": used, "obs": obs, "xmap": {c: i for i, c in enumerate(used)}, "xlabel": "condition",
                        "title": obs.replace("__", " ")})
    return out


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
    for ax, panel in zip(axes, panels):
        ax.set_facecolor(SURFACE)
        obs, xmap = panel["obs"], panel["xmap"]
        series = panel.get("series")                      # one ELISA series of an observable (None: all of them)
        sel = meas["simulationConditionId"].isin(panel["conds"]) & ((meas["observableId"] == obs) if obs else True)
        if series is not None:
            sel &= meas["observableParameters"] == series
        d = meas[sel]
        xof = (lambda t: t["time"]) if xmap is None else (lambda t: t["simulationConditionId"].map(xmap))
        for key in ["published"] + order[::-1]:           # worst route first, widest: the best is the thin line on top
            c = curves[key]
            cs = c[c["simulationConditionId"].isin(panel["conds"]) & ((c["observableId"] == obs) if obs else True)]
            if series is not None:
                cs = cs[cs["observableParameters"] == series]
            for _o, g in (cs.groupby(["observableId", "simulationConditionId"]) if xmap is None else [(None, cs)]):   # one line per condition
                g = g.assign(_x=xof(g)).sort_values("_x")
                mark = dict(marker="o", markersize=7, markerfacecolor="none") if len(g) == 1 or xmap is not None else {}   # points, not a line
                categorical = xmap is not None and panel["xlabel"] == "condition"       # separate conditions: points, no line
                if key == "published":
                    ax.plot(g["_x"], g["simulation"], color="#8a8985", lw=1.3, ls="none" if categorical else (0, (4, 3)), zorder=2, **mark)
                else:
                    rank = order.index(key)
                    ax.plot(g["_x"], g["simulation"], color=colors[key], zorder=3, ls="none" if categorical else "-",
                            lw=1.4 + 1.8 * rank / max(len(order) - 1, 1), **mark)
        ax.scatter(xof(d), d["measurement"], s=16, color=INK, zorder=5, edgecolor=SURFACE, linewidth=0.6)
        if xmap is not None and panel["xlabel"] == "condition":
            ax.set_xticks(range(len(panel["conds"])), panel.get("ticks") or [c.split("__")[-1] for c in panel["conds"]],
                          fontsize=6.5)
            ax.set_xlim(-0.5, len(panel["conds"]) - 0.5)
        elif xmap is not None and (lambda pos: len(pos) > 1 and max(pos) / min(pos) > 100)([v for v in xmap.values() if v > 0]):
            if min(xmap.values()) > 0:
                ax.set_xscale("log")
            else:                                         # a dose of zero is there too (Schwen's IR sum)
                ax.set_xscale("symlog", linthresh=min(v for v in xmap.values() if v > 0))
        if (problem.observables["observableTransformation"].reindex([obs]).iloc[0] if obs and "observableTransformation" in problem.observables
                else "lin") in ("log", "log10"):
            ax.set_yscale("log")                          # the observable is fitted on a log scale, so it is drawn on one
        ax.set_title(panel["title"], loc="left", fontsize=9.5, color=INK)
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
    for ax, panel in zip(axes, panels):
        if panel.get("xlabel") != "condition":
            ax.set_xlabel(panel.get("xlabel", "time"), color=MUTED, fontsize=8, labelpad=2)
    fig.suptitle(f"{args.problem}: final fit of every route, re-simulated with PEtab's reference simulator",
                 x=0.01, ha="left", fontsize=12, color=INK)
    fig.text(0.01, 0.005, f"Lower NLL is better; routes within {spread:.2f} NLL of each other overlay. One start per curve, "
             "so routes may end in different optima." if spread < 5 else
             "Lower NLL is better. One start per curve: routes end in different optima" + (
                 ", and the published parameters were not fitted with this objective's free noise parameters."
                 if noise_names(problem) - {"sigma"} else ", and the published parameters are PEtab's nominal values, not a fit."), fontsize=8, color=MUTED, va="bottom")
    fig.tight_layout(rect=[0, 0.035 + (0.025 + 0.018 * math.ceil((len(order) + len(failed) + 2) / 2) if long_legend else 0), 1, 0.97])
    fig.savefig(out / "fits.png", dpi=150, facecolor=SURFACE)
    print("wrote", out / "fits.png")
    for k in order:
        print(f"  {nll[k]:+10.4f}  {labels[k]}")
    print(f"  {nll['published']:+10.4f}  {labels['published']}")


if __name__ == "__main__":
    main()

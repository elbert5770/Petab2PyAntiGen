"""A PEtab problem straight through AMICI, without PyAntiGen: the baseline PyAntiGen's AMICI path is measured against.

    python benchmarks/amici_direct.py --problem Zhao_QuantBiol2020 --spec OPTIMIZATION_Zhao

Imports the untouched PEtab problem with AMICI's PEtab importer (compiled once), then runs L-BFGS-B on the log10 parameters
with AMICI's likelihood and sensitivities (``simulate_petab``), from the same perturbed start as
``benchmarks/fit_problem.py`` (seed 1, x * 10**U(-0.3, 0.3) in the PyAntiGen spec's parameter order). Every parameter PEtab
estimates is free, including its noise parameters (``sigma``, ``sd_*``), which start at their maximum-likelihood value for
the starting parameters: started at their nominal values while the residuals are far larger, the joint (theta, sigma)
problem is badly conditioned and L-BFGS-B stops early at a poor point. PyAntiGen's profiled sigma has no such problem.

Run in WSL/Linux with AMICI and the ``petab`` package (a venv layered on AMICI's, plus ``pip install petab``; see the
README, "AMICI (WSL)"):

    PATH=/opt/venv/bin:$PATH python benchmarks/amici_direct.py --problem Zheng_PNAS2012 --spec OPTIMIZATION_Zheng

Paths: PYANTIGEN_DIR, PETAB_BENCHMARK_DIR, AMICI_OUT override the defaults (sibling checkouts, ~/amici_direct).
Writes benchmarks/<problem>/fits/direct_amici_lbfgsb.json (``..._seed<n>.json`` for a start other than seed 1).
"""
import argparse
import json
import logging
import math
import os
import platform
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [os.environ.get("PYANTIGEN_DIR", str(ROOT.parent / "PyAntiGen")), str(ROOT), str(HERE)]

import numpy as np
from scipy.optimize import minimize
import petab.v1 as petab
import amici
from amici.importers.petab.v1 import import_petab_problem
from amici.sim.sundials.petab.v1 import simulate_petab, LLH, SLLH

logging.getLogger("amici").setLevel(logging.ERROR)
for name in list(logging.root.manager.loggerDict):
    if name.startswith("amici"):
        logging.getLogger(name).setLevel(logging.ERROR)

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--problem", default="Zheng_PNAS2012")
ap.add_argument("--spec", required=True, help="the PyAntiGen Optimization whose parameters and order define the start")
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--yaml", default=None, help="a variant of the PEtab problem (problem.yaml) to fit instead of the benchmark's")
ap.add_argument("--out", default=None, help="folder for the result (default benchmarks/<problem>/fits)")
ap.add_argument("--max-evals", type=int, default=6000, help="evaluation budget (L-BFGS-B maxfun)")
ap.add_argument("--rtol", type=float, default=None, help="solver relative tolerance (default: AMICI's)")
ap.add_argument("--atol", type=float, default=None, help="solver absolute tolerance (default: AMICI's, 1e-16)")
ap.add_argument("--retry", type=lambda v: [tuple(float(x) for x in pair.split(":")) for pair in v.split(",")], default=[],
                help="looser rtol:atol pairs to retry a failed integration with, e.g. 1e-6:1e-8,1e-5:1e-7")
args = ap.parse_args()

# --- the PyAntiGen spec: parameter names, order and bounds, so the start is the same point
project = ROOT / "problems" / args.problem
os.chdir(project)
sys.path.insert(0, str(project))
import AntiGen_paths                                                     # noqa: F401,E402
import Model_run                                                         # noqa: E402
from fit_problem import jittered_start                                  # noqa: E402

spec = Model_run.get_OPTIMIZATION(args.spec)
names = list(spec.param_names)
start = jittered_start(spec.x0, spec.bounds if isinstance(spec.bounds, dict) else None, args.seed)

# --- the PEtab problem, through AMICI
bench = Path(os.environ.get("PETAB_BENCHMARK_DIR", ROOT.parent / "petab-benchmark-problems"))
YAML = Path(args.yaml) if args.yaml else bench / "problems" / args.problem / "v1" / "problem.yaml"
prob = petab.Problem.from_yaml(str(YAML))
t = time.time()
model = import_petab_problem(prob, output_dir=os.environ.get("AMICI_OUT", os.path.expanduser(f"~/amici_direct/{args.problem}")),
                             model_name=re.sub(r"\W", "_", args.problem) + "_direct")
print(f"RESULT import/compile: {time.time() - t:.1f} s")
solver = model.create_solver()
solver.set_sensitivity_order(amici.sim.sundials.SensitivityOrder.first)
solver.set_max_steps(2000000)           # AMICI's default step budget is small for a stiff network (PyAntiGen's default is 2e6 too)
if args.rtol:
    solver.set_relative_tolerance(args.rtol)
if args.atol:       # Smith: the sensitivity integration fails ("too much work") at AMICI's default atol of 1e-16
    solver.set_absolute_tolerance(args.atol)

pt = prob.parameter_df
assert set(pt.loc[pt["estimate"] == 1, "parameterScale"]) <= {"log10", "lin"}, "this script fits log10 and linear parameters"
nominal = dict(zip(pt.index, pt["nominalValue"]))
estimated = list(pt.index[pt["estimate"] == 1])
free = names + [p for p in estimated if p not in names]                  # then the noise parameters
is_log = {p: pt.loc[p, "parameterScale"] == "log10" for p in free}       # a "lin" parameter (Schwen: fragments, scaleElisa_*) is fitted as it is
to_x = lambda p, v: math.log10(v) if is_log[p] else float(v)           # noqa: E731
from_x = lambda p, x: 10 ** x if is_log[p] else float(x)               # noqa: E731
lo = {p: to_x(p, float(pt.loc[p, "lowerBound"])) for p in free}
hi = {p: to_x(p, float(pt.loc[p, "upperBound"])) for p in free}
N = len(prob.measurement_df)

# objective priors (Schwen: parameterScaleNormal on six rate constants): AMICI's likelihood does not include them, PEtab's objective does.
# nll += 0.5 ((x - mean) / sd)^2 for x the parameter on its scale; d/d(log10 p) = (x - mean) / sd^2 for a log10 parameter.
PRIORS = {}
if "objectivePriorType" in pt:
    for p in free:
        kind = pt.loc[p, "objectivePriorType"]
        if isinstance(kind, str) and kind:
            assert kind == "parameterScaleNormal", f"objective prior {kind} on {p} is not implemented"
            PRIORS[p] = tuple(float(v) for v in str(pt.loc[p, "objectivePriorParameters"]).split(";"))


def prior_value_and_gradient(x_of):
    """(nll, {parameter: d nll / d x}) of the priors at x_of[p] (the parameter on its scale)."""
    nll, grad = 0.0, {}
    for p, (mean, sd) in PRIORS.items():
        nll += 0.5 * ((x_of[p] - mean) / sd) ** 2
        grad[p] = (x_of[p] - mean) / sd ** 2
    return nll, grad


# a measurement on a log scale carries the Jacobian of its transformation in PEtab's density: a constant of the data
_trafo = prob.observable_df.get("observableTransformation")
JAC = 0.0
if _trafo is not None:
    for _o, _m in zip(prob.measurement_df["observableId"], prob.measurement_df["measurement"]):
        _k = str(_trafo.get(_o, "lin"))
        JAC += math.log(_m) if _k == "log" else math.log(_m * math.log(10)) if _k == "log10" else 0.0
# PyAntiGen prints the NLL without the Gaussian constants: N/2 ln(2 pi) for every point, and the N_est/2 that profiling a sigma
# leaves for the points of an ESTIMATED noise parameter (a declared sigma, estimate = 0, leaves n ln(sigma) + SSE/(2 sigma^2))
def _is_number(x):
    try:
        float(x)
        return True
    except (TypeError, ValueError):
        return False


if "noiseParameters" in prob.measurement_df:
    N_EST = int(prob.measurement_df["noiseParameters"].isin(estimated).sum())
else:                                       # sigma is in the observable table: a constant (Smith: 1.0), or a free parameter
    N_EST = 0 if all(_is_number(f) for f in prob.observable_df["noiseFormula"]) else N
CONST = 0.5 * N * math.log(2 * math.pi) + 0.5 * N_EST + JAC          # (the priors' own ln(sd sqrt(2 pi)) is not in the objective here)


def evaluate(pars):
    r = simulate_petab(prob, amici_model=model, solver=solver, problem_parameters=pars, scaled_parameters=False)
    if r[SLLH] is None and args.retry:        # an integration failed (Smith: a Heaviside root, knife-edge in the parameters):
        keep = (solver.get_relative_tolerance(), solver.get_absolute_tolerance())       # retry at looser tolerances, as PyAntiGen's
        for rtol, atol in args.retry:                                                   # analytic_gradient "retry" does
            solver.set_relative_tolerance(rtol)
            solver.set_absolute_tolerance(atol)
            r = simulate_petab(prob, amici_model=model, solver=solver, problem_parameters=pars, scaled_parameters=False)
            if r[SLLH] is not None:
                break
        solver.set_relative_tolerance(keep[0])
        solver.set_absolute_tolerance(keep[1])
        calls["retried"] += 1
    return r


# noise parameters start at their ML value: one per noiseParameters entry of the measurement table
# (llh = -n/2 ln(2 pi sigma^2) - SSE/(2 sigma^2), so SSE follows from the llh at sigma = 1)
start_pars = {**nominal, **start}
noise = [p for p in free if p not in names]
groups = {p: (prob.measurement_df["noiseParameters"] == p) for p in noise}
for g in noise:
    start_pars[g] = 1.0
r1 = evaluate(start_pars)
sse_total = 2.0 * (-r1[LLH] - 0.5 * N * math.log(2 * math.pi) - JAC)         # all datasets at sigma = 1
if not noise:                                                            # every sigma is declared (Smith): nothing to start
    pass
elif len(noise) == 1:
    start_pars[noise[0]] = math.sqrt(sse_total / N)
else:                                                                    # per-dataset: use the PEtab reference's SSE per group
    from petab2pyantigen.petab_reader import load_petab
    from petab2pyantigen.reference import _transform, simulate_petab as reference_simulate
    sims = reference_simulate(load_petab(str(YAML)), {k: v for k, v in start_pars.items() if k in names})
    for g in noise:
        rows = sims[sims["noiseParameters"] == g]
        kinds = rows["observableId"].map(lambda o: str(_trafo.get(o, "lin")) if _trafo is not None else "lin")
        res = [_transform(m, k) - _transform(y, k) for m, y, k in zip(rows["measurement"], rows["simulation"], kinds)]
        start_pars[g] = min(max(math.sqrt(float(np.sum(np.square(res))) / len(rows)), float(pt.loc[g, "lowerBound"])),
                            float(pt.loc[g, "upperBound"]))             # PEtab bounds the noise parameter
print("RESULT start: noise parameters at their ML values:", {g: round(start_pars[g], 4) for g in noise[:3]}, "...")

calls = {"n": 0, "t": 0.0, "retried": 0}


def fun(x):
    pars = dict(nominal)
    pars.update({n: from_x(n, v) for n, v in zip(free, x)})
    t = time.time()
    r = evaluate(pars)
    calls["t"] += time.time() - t
    calls["n"] += 1
    if r[SLLH] is None:                      # an integration failed: say which condition (the likelihood is nan then)
        from amici.sim.sundials.petab.v1 import RDATAS
        bad = [(rd.id, rd.status) for rd in r[RDATAS] if rd.status != 0][:3]
        raise RuntimeError(f"AMICI could not integrate (status per condition: {bad}); try --rtol/--atol, a looser pair")
    g = np.array([-r[SLLH][n] * (pars[n] * math.log(10) if is_log[n] else 1.0) for n in free])     # d(-llh)/d x, x = the fitted scale
    nll = -r[LLH]
    if PRIORS:
        pn, pg = prior_value_and_gradient({p: to_x(p, pars[p]) for p in PRIORS})
        nll += pn
        g = g + np.array([pg.get(n, 0.0) for n in free])
    if calls["n"] % 25 == 0:
        print(f"RESULT progress: {calls['n']} calls, -llh - constant = {nll - CONST:.4f}, {time.time() - t0:.0f} s", flush=True)
    return nll, g


x0 = np.array([to_x(n, start_pars[n]) for n in free])
t0 = time.time()
res = minimize(fun, x0, jac=True, method="L-BFGS-B", bounds=[(lo[p], hi[p]) for p in free],
               options={"maxiter": 2000, "maxfun": args.max_evals})
wall = time.time() - t0
print(f"RESULT success={res.success} msg={res.message}")
print(f"RESULT integration failures retried at a looser tolerance: {calls['retried']} of {calls['n']} calls")
print(f"RESULT nit={res.nit} nfev={res.nfev} wall={wall:.1f}s  in simulate_petab={calls['t']:.1f}s ({1000 * calls['t'] / calls['n']:.0f} ms/call)")
print(f"RESULT -llh={res.fun:.4f}  NLL without constant={res.fun - CONST:.4f}")

out = Path(args.out) if args.out else HERE / args.problem / "fits"
out.mkdir(parents=True, exist_ok=True)
label = "direct_amici_lbfgsb" + ("" if args.seed == 1 else f"_seed{args.seed}")
(out / f"{label}.json").write_text(json.dumps({
    "problem": args.problem, "label": label, "route": "PEtab -> AMICI directly / L-BFGS-B", "backend": "amici",
    "method": "L-BFGS-B", "sigma": "estimated jointly (PEtab parameters)", "exact_derivatives": True, "seed": args.seed,
    "params": {n: float(from_x(n, v)) for n, v in zip(free, res.x)},
    "objective": float(res.fun - CONST), "success": bool(res.success), "message": str(res.message),
    "nfev": int(res.nfev), "nit": int(res.nit), "wall_s": round(wall, 1),
    "platform": platform.platform(), "python": platform.python_version(),
    "note": "objective = -llh + priors - N/2 ln(2 pi) - N_est/2 (N_est: points of estimated noise parameters)"
            + ("; PEtab's objective priors added (AMICI's llh has none)" if PRIORS else "")
            + (f"; solver rtol {args.rtol}, atol {args.atol}, retry {args.retry}" if args.rtol or args.atol else ""),
}, indent=1), encoding="utf-8")

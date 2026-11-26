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
Writes benchmarks/<problem>/fits/direct_amici_lbfgsb.json.
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
YAML = bench / "problems" / args.problem / "v1" / "problem.yaml"
prob = petab.Problem.from_yaml(str(YAML))
t = time.time()
model = import_petab_problem(prob, output_dir=os.environ.get("AMICI_OUT", os.path.expanduser(f"~/amici_direct/{args.problem}")),
                             model_name=re.sub(r"\W", "_", args.problem) + "_direct")
print(f"RESULT import/compile: {time.time() - t:.1f} s")
solver = model.create_solver()
solver.set_sensitivity_order(amici.sim.sundials.SensitivityOrder.first)

pt = prob.parameter_df
assert (pt.loc[pt["estimate"] == 1, "parameterScale"] == "log10").all(), "this script fits log10-scaled parameters"
nominal = dict(zip(pt.index, pt["nominalValue"]))
estimated = list(pt.index[pt["estimate"] == 1])
free = names + [p for p in estimated if p not in names]                  # then the noise parameters
lo = {p: math.log10(float(pt.loc[p, "lowerBound"])) for p in free}
hi = {p: math.log10(float(pt.loc[p, "upperBound"])) for p in free}
N = len(prob.measurement_df)
CONST = 0.5 * N * (1 + math.log(2 * math.pi))                            # PyAntiGen prints NLL without this constant


def evaluate(pars):
    return simulate_petab(prob, amici_model=model, solver=solver, problem_parameters=pars, scaled_parameters=False)


# noise parameters start at their ML value: one per noiseParameters entry of the measurement table
# (llh = -n/2 ln(2 pi sigma^2) - SSE/(2 sigma^2), so SSE follows from the llh at sigma = 1)
start_pars = {**nominal, **start}
noise = [p for p in free if p not in names]
groups = {p: (prob.measurement_df["noiseParameters"] == p) for p in noise}
for g in noise:
    start_pars[g] = 1.0
r1 = evaluate(start_pars)
sse_total = 2.0 * (-r1[LLH] - 0.5 * N * math.log(2 * math.pi))               # all datasets at sigma = 1
if len(noise) == 1:
    start_pars[noise[0]] = math.sqrt(sse_total / N)
else:                                                                    # per-dataset: use the PEtab reference's SSE per group
    from petab2pyantigen.petab_reader import load_petab
    from petab2pyantigen.reference import simulate_petab as reference_simulate
    sims = reference_simulate(load_petab(str(YAML)), {k: v for k, v in start_pars.items() if k in names})
    for g in noise:
        rows = sims[sims["noiseParameters"] == g]
        start_pars[g] = math.sqrt(float(((rows["measurement"] - rows["simulation"]) ** 2).sum()) / len(rows))
print("RESULT start: noise parameters at their ML values:", {g: round(start_pars[g], 4) for g in noise[:3]}, "...")

calls = {"n": 0, "t": 0.0}


def fun(x):
    pars = dict(nominal)
    pars.update({n: 10 ** v for n, v in zip(free, x)})
    t = time.time()
    r = evaluate(pars)
    calls["t"] += time.time() - t
    calls["n"] += 1
    g = np.array([-r[SLLH][n] * pars[n] * math.log(10) for n in free])         # d(-llh)/d log10 p
    return -r[LLH], g


x0 = np.array([math.log10(start_pars[n]) for n in free])
t0 = time.time()
res = minimize(fun, x0, jac=True, method="L-BFGS-B", bounds=[(lo[p], hi[p]) for p in free],
               options={"maxiter": 2000, "maxfun": 6000})
wall = time.time() - t0
print(f"RESULT success={res.success} msg={res.message}")
print(f"RESULT nit={res.nit} nfev={res.nfev} wall={wall:.1f}s  in simulate_petab={calls['t']:.1f}s ({1000 * calls['t'] / calls['n']:.0f} ms/call)")
print(f"RESULT -llh={res.fun:.4f}  NLL without constant={res.fun - CONST:.4f}")

out = HERE / args.problem / "fits"
out.mkdir(parents=True, exist_ok=True)
(out / "direct_amici_lbfgsb.json").write_text(json.dumps({
    "problem": args.problem, "label": "direct_amici_lbfgsb", "route": "PEtab -> AMICI directly / L-BFGS-B", "backend": "amici",
    "method": "L-BFGS-B", "sigma": "estimated jointly (PEtab parameters)", "exact_derivatives": True, "seed": args.seed,
    "params": {n: 10 ** float(v) for n, v in zip(free, res.x)},
    "objective": float(res.fun - CONST), "success": bool(res.success), "message": str(res.message),
    "nfev": int(res.nfev), "nit": int(res.nit), "wall_s": round(wall, 1),
    "platform": platform.platform(), "python": platform.python_version(),
    "note": "objective = -llh - n/2 (1 + ln 2 pi)",
}, indent=1), encoding="utf-8")

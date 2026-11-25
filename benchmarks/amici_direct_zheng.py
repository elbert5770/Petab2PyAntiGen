"""A PEtab problem straight through AMICI, without PyAntiGen: the baseline PyAntiGen's AMICI path is measured against.

Imports the untouched PEtab problem with AMICI's PEtab importer (compiled once, ~15 s for Zheng), then runs
L-BFGS-B on the log10 parameters with AMICI's likelihood and sensitivities (``simulate_petab``), from the same
perturbed start as ``Model_run.py --optimize ZhengAMICI`` (seed 1, x * 10**U(-0.3, 0.3), PyAntiGen's parameter order).
PEtab's explicit ``sigma`` starts at its ML value for the starting parameters; started at its nominal 0.0023 while the
residuals are ~0.1, the joint (theta, sigma) problem is badly conditioned and L-BFGS-B stops after 7 iterations at sigma = 0.11.

Run in WSL/Linux with AMICI and the ``petab`` package (this repository's notes: README, "AMICI (WSL)"):

    PATH=/opt/venv/bin:$PATH python benchmarks/amici_direct_zheng.py

``petab`` is not in AMICI's own environment: use a venv layered on it (a .pth file pointing at its site-packages) with
``pip install petab``. Paths: PYANTIGEN_DIR, PETAB_YAML, AMICI_OUT override the defaults (sibling checkouts, ~/amici_direct).

Measured on this machine (Zheng_PNAS2012, 46 parameters): 14.6 s import, 48 ms per likelihood, 113 ms per likelihood +
gradient, fit converged in 474 evaluations / 58 s to -llh = -278.3286 (PEtab nominal -278.3335).
"""

import os, sys, time, logging, math
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))          # this repository
sys.path[:0] = [os.environ.get("PYANTIGEN_DIR", os.path.join(os.path.dirname(HERE), "PyAntiGen")), HERE]
import numpy as np
from scipy.optimize import minimize
import petab.v1 as petab
import amici
from amici.importers.petab.v1 import import_petab_problem
from amici.sim.sundials.petab.v1 import simulate_petab, LLH, SLLH
from antimony_modules.Histone.methylation_Zheng_PNAS2012 import rate_constants

logging.getLogger("amici").setLevel(logging.ERROR)
for name in list(logging.root.manager.loggerDict):
    if name.startswith("amici"):
        logging.getLogger(name).setLevel(logging.ERROR)

YAML = os.environ.get("PETAB_YAML", os.path.join(os.path.dirname(HERE), "petab-benchmark-problems", "problems",
                                                  "Zheng_PNAS2012", "v1", "problem.yaml"))
prob = petab.Problem.from_yaml(YAML)
model = import_petab_problem(prob, output_dir=os.environ.get("AMICI_OUT", os.path.expanduser("~/amici_direct/Zheng")), model_name="Zheng_direct")
solver = model.create_solver()
solver.set_sensitivity_order(amici.sim.sundials.SensitivityOrder.first)

nominal = dict(zip(prob.parameter_df.index, prob.parameter_df["nominalValue"]))
names = ["inflowp"] + rate_constants()                       # PyAntiGen's order: the jitter must hit the same parameters
rng = np.random.default_rng(1)
start = dict(nominal)
for n in names:
    start[n] = float(np.clip(nominal[n] * 10 ** rng.uniform(-0.3, 0.3), 1e-5, 1e3))
free = names + ["sigma"]                                     # sigma stays at its nominal start
CONST = 0.5 * 60 * (1 + math.log(2 * math.pi))               # PyAntiGen prints NLL without this constant

calls = {"n": 0, "t": 0.0}
best = {"f": np.inf, "n": 0, "t": 0.0}
t0 = time.time()


def fun(x):
    pars = dict(nominal)
    pars.update({n: 10 ** v for n, v in zip(free, x)})
    t = time.time()
    r = simulate_petab(prob, amici_model=model, solver=solver, problem_parameters=pars, scaled_parameters=False)
    calls["t"] += time.time() - t
    calls["n"] += 1
    f = -r[LLH]
    g = np.array([-r[SLLH][n] * pars[n] * math.log(10) for n in free])      # d(-llh)/d log10 p
    if f < best["f"]:
        best.update(f=f, n=calls["n"], t=time.time() - t0)
    return f, g


# start sigma at its ML value for the starting parameters: llh(sigma=1) = -30 log(2 pi) - SSE/2
_p = dict(start); _p["sigma"] = 1.0
_r = simulate_petab(prob, amici_model=model, solver=solver, problem_parameters=_p, scaled_parameters=False)
_sse = 2.0 * (-_r[LLH] - 30 * math.log(2 * math.pi))
start["sigma"] = math.sqrt(_sse / 60)
print(f"RESULT start: SSE={_sse:.4g}  sigma_hat={start['sigma']:.4g}")
x0 = np.array([math.log10(start[n]) for n in free])
res = minimize(fun, x0, jac=True, method="L-BFGS-B", bounds=[(-5, 3)] * len(free),
               options={"maxiter": 1000, "maxfun": 3000})
wall = time.time() - t0
print(f"RESULT success={res.success} msg={res.message}")
print(f"RESULT nit={res.nit} nfev={res.nfev} wall={wall:.1f}s  in simulate_petab={calls['t']:.1f}s ({1000*calls['t']/calls['n']:.0f} ms/call)")
print(f"RESULT -llh={res.fun:.4f}  NLL without constant={res.fun - CONST:.4f}  (PEtab nominal: -278.3335 / -363.4700)")
print(f"RESULT sigma={10**res.x[-1]:.5f}")

# the same record PyAntiGen's fit_zheng.py writes, for compare_fits.py
import json, platform
out = os.path.join(HERE, "benchmarks", "fits")
os.makedirs(out, exist_ok=True)
with open(os.path.join(out, "direct_amici_lbfgsb.json"), "w", encoding="utf-8") as f:
    json.dump({
        "label": "direct_amici_lbfgsb", "route": "PEtab -> AMICI directly / L-BFGS-B", "backend": "amici",
        "method": "L-BFGS-B", "sigma": "estimated jointly (PEtab parameter)", "exact_derivatives": True, "seed": 1,
        "params": {n: 10 ** float(v) for n, v in zip(free, res.x)},
        "objective": float(res.fun - CONST), "success": bool(res.success), "message": str(res.message),
        "nfev": int(res.nfev), "nit": int(res.nit), "wall_s": round(wall, 1),
        "platform": platform.platform(), "python": platform.python_version(), "note": "objective = -llh - n/2 (1 + ln 2 pi)",
    }, f, indent=1)

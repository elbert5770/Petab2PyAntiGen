"""Check a PyAntiGen project against the PEtab problem it replaces.

Both sides simulate the same experiment -- PEtab through ``reference.simulate_petab``
(SBML + tables, RoadRunner steady state), PyAntiGen through its own engine
(generated Antimony, simulation blocks, events) -- and the observables are compared
at every measurement time, at the nominal parameters and at random ones.

    python -m petab2pyantigen.verify --petab <PEtab dir> --project problems/<Name> \
           --experiment <key in EXPERIMENT_dict> [--random 3]
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from .petab_reader import load_petab
from .reference import log_likelihood, nominal_values, simulate_petab

_DRIVER = r'''
import json, os, pickle, sys
sys.path.insert(0, os.getcwd())
import AntiGen_paths
from Model_run import EXPERIMENT_dict, MODEL_NAME, REPO_ROOT
from pyantigen.engine.Model_simulate import setup_simulation

out, key, overrides = sys.argv[1], sys.argv[2], json.loads(sys.argv[3])
captured = {}
def capture(paths, results_dict):
    for label, item in results_dict.items():
        res = item["results"]
        captured[label] = {"results": {c: res[:, i].copy() for i, c in enumerate(res.colnames)},
                           "data": item["data"]}
exp = dict(EXPERIMENT_dict[key]); exp["plot"] = capture
settings = {"run_steady_state_first": False, "Verbose": False, "save_SBML?": False,
            "MODEL_NAME": MODEL_NAME, "REPO_ROOT": REPO_ROOT}
# Apply the overrides the way a fit applies a parameter vector: after the replicate's Update_parameters, set what is a model
# parameter, then the replicate's parameter_hooks and Update_opt_parameters. (A fitted name may only mean something to a
# hook: R_Stage_I_Wuhan sets R in the Wuhan, Stage I replicate and does nothing in the others.)
if overrides:
    from pyantigen.engine.Optimize import set_parameters_from_dict
    for rep in exp["EXPERIMENT"].replicates.values():
        def _update(r, replicate, _first=rep["Update_parameters"]):
            _first(r, replicate)
            set_parameters_from_dict(r, overrides)
            for hook in replicate.get("parameter_hooks", []):
                hook(r, overrides)
            if replicate.get("Update_opt_parameters") is not None:
                replicate["Update_opt_parameters"](r, replicate, overrides)
        rep["Update_parameters"] = _update
setup_simulation(settings, exp)
pickle.dump(captured, open(out, "wb"))
'''


def run_pyantigen(project_dir: str | Path, experiment: str, parameters: dict[str, float] | None = None) -> dict:
    """Simulate a project with its own engine; returns {replicate label: {"results": {col: array}, "data": df}}."""
    project_dir = Path(project_dir).resolve()
    with tempfile.TemporaryDirectory() as tmp:
        driver, out = Path(tmp) / "driver.py", Path(tmp) / "out.pkl"
        driver.write_text(_DRIVER, encoding="utf-8")
        proc = subprocess.run([sys.executable, "-I", str(driver), str(out), experiment, json.dumps(parameters or {})],
                              cwd=project_dir, capture_output=True, text=True,
                              env={**os.environ, "PYTHONPATH": ""})
        if proc.returncode != 0:
            raise RuntimeError(f"PyAntiGen simulation failed:\n{proc.stdout[-2000:]}\n{proc.stderr[-3000:]}")
        return pickle.loads(out.read_bytes())


def compare(petab_dir, project_dir, experiment, parameters=None, observable_name=None, strip_prefix="",
            replicate_of=None) -> dict:
    """Max abs/rel difference between PEtab's and PyAntiGen's observables at every measurement.

    ``observable_name`` maps a PEtab observable id to the PyAntiGen result column (``{id: column}`` or a callable;
    default: the same name). A PEtab condition is a PyAntiGen replicate: ``replicate_of(condition_id)``, by default the
    condition id without ``strip_prefix`` (or, with one replicate, that replicate).
    """
    problem = load_petab(petab_dir)
    ref = simulate_petab(problem, parameters)
    pa = run_pyantigen(project_dir, experiment, parameters)
    results = {label: pd.DataFrame(item["results"]) for label, item in pa.items()}
    if replicate_of is None:
        if len(pa) == 1:
            replicate_of = lambda _c: next(iter(pa))                                          # noqa: E731
        else:
            replicate_of = lambda c: c[len(strip_prefix):] if c.startswith(strip_prefix) else c   # noqa: E731
    column = (observable_name if callable(observable_name)
              else (lambda o, m=observable_name or {}: m.get(o, o)))
    # the time at which each replicate's data start, against the PEtab time of its condition (a pre-equilibration
    # block shifts one by its length)
    offset = {}
    for cond in problem.measurements["simulationConditionId"].unique():
        label = replicate_of(cond)
        petab_t0 = problem.measurements.loc[problem.measurements["simulationConditionId"] == cond, "time"].min()
        offset[cond] = float(pa[label]["data"][label]["time"].min() - petab_t0)
    got = []
    for _, m in ref.iterrows():
        label, res = replicate_of(m["simulationConditionId"]), None
        res = results[label]
        want = m["time"] + offset[m["simulationConditionId"]]
        i = int(np.argmin(np.abs(res["time"].values - want)))
        if abs(res["time"].values[i] - want) > 1e-6 * max(1.0, abs(want)):
            raise RuntimeError(f"no PyAntiGen output at t = {want} in {label}")
        got.append(res[column(m["observableId"])].values[i])
    ref = ref.assign(pyantigen=got)
    diff = (ref["pyantigen"] - ref["simulation"]).abs()
    return {"n": len(ref), "max_abs": float(diff.max()),
            "max_rel": float((diff / ref["simulation"].abs().clip(lower=1e-300)).max()),
            "llh_petab": log_likelihood(problem, ref.assign(simulation=ref["simulation"]), parameters),
            "llh_pyantigen": log_likelihood(problem, ref.assign(simulation=ref["pyantigen"]), parameters),
            "table": ref}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--petab", required=True, help="PEtab v1 directory (contains problem.yaml)")
    ap.add_argument("--project", required=True, help="PyAntiGen project folder (problems/<Name>)")
    ap.add_argument("--experiment", required=True, help="key of EXPERIMENT_dict in the project's Model_run.py")
    ap.add_argument("--random", type=int, default=3, help="number of random parameter sets (log-uniform x0.3..x3)")
    ap.add_argument("--tol", type=float, default=1e-6, help="max abs difference accepted")
    ap.add_argument("--rel-tol", type=float, default=None,
                    help="accept a relative difference instead (for observables that span orders of magnitude)")
    ap.add_argument("--strip-prefix", default="", help="PyAntiGen replicate label = PEtab condition id without this prefix")
    ap.add_argument("--observable", action="append", default=[], metavar="PETAB_ID=COLUMN",
                    help="PEtab observable -> PyAntiGen result column (repeatable; default: the same name)")
    args = ap.parse_args(argv)
    obs_map = dict(item.split("=", 1) for item in args.observable)

    problem = load_petab(args.petab)
    nominal = nominal_values(problem)
    est = [p for p, e in problem.parameters["estimate"].items() if e == 1]
    rng = np.random.default_rng(0)
    lo, hi = problem.parameters["lowerBound"].astype(float), problem.parameters["upperBound"].astype(float)
    cases = [("nominal", None)] + [
        (f"random {i + 1}", {p: float(np.clip(nominal[p] * 10 ** rng.uniform(-0.5, 0.5), lo[p], hi[p])) for p in est
                             if p in set(problem.parameters.index)})
        for i in range(args.random)]
    ok = True
    for name, pars in cases:
        r = compare(args.petab, args.project, args.experiment, pars, observable_name=obs_map,
                    strip_prefix=args.strip_prefix)
        good = (r["max_rel"] <= args.rel_tol) if args.rel_tol is not None else (r["max_abs"] <= args.tol)
        ok &= good
        print(f"{name:10s} n={r['n']}  max|diff|={r['max_abs']:.3e}  max rel={r['max_rel']:.3e}  "
              f"llh PEtab={r['llh_petab']:.6f}  llh PyAntiGen={r['llh_pyantigen']:.6f}  {'OK' if good else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())

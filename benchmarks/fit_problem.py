"""Fit a converted PEtab problem with one PyAntiGen route and save the result in the format compare_fits.py reads.

    python benchmarks/fit_problem.py --problem Zheng_PNAS2012 --label pyantigen_amici_lbfgsb --spec OPTIMIZATION_Zheng_amici

``--spec`` names an Optimization in problems/<problem>/Modules/Optimizer_settings.py. Every route starts from the same
perturbed point (seed 1: each parameter times 10**U(-0.3, 0.3), kept inside its bounds, in the spec's parameter order),
the start benchmarks/amici_direct.py uses too. Writes benchmarks/<problem>/fits/<label>.json.
"""
import argparse
import json
import os
import platform
import shutil
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def jittered_start(x0: dict, bounds: dict, seed: int) -> dict:
    """x * 10**U(-0.3, 0.3), clipped to each parameter's bounds; seed 0 keeps the nominal point."""
    if not seed:
        return dict(x0)
    rng = np.random.default_rng(seed)
    out = {}
    for k, v in x0.items():
        lo, hi = bounds[k] if bounds and k in bounds else (-np.inf, np.inf)
        out[k] = float(np.clip(v * 10 ** rng.uniform(-0.3, 0.3), lo, hi))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--problem", default="Zheng_PNAS2012", help="folder under problems/")
    ap.add_argument("--label", required=True, help="name of the result file (benchmarks/<problem>/fits/<label>.json)")
    ap.add_argument("--spec", required=True, help="Optimization in Modules/Optimizer_settings.py")
    ap.add_argument("--seed", type=int, default=1, help="start perturbation (0: PEtab's nominal point)")
    ap.add_argument("--max-evals", type=int, default=None, help="evaluation budget (default: the spec's)")
    ap.add_argument("--backend", default=None, choices=["roadrunner", "amici", "auto"])
    ap.add_argument("--no-exact", action="store_true", help="finite differences of the whole loss instead of the chain rule")
    ap.add_argument("--observables", choices=["column", "callable"], default="column",
                    help="column: bare result columns; callable: the same numbers from a function of the table, which the "
                         "engine cannot take as a bare column (as a difference of simulations is). Zheng only.")
    ap.add_argument("--note", default="")
    args = ap.parse_args(argv)

    project = ROOT / "problems" / args.problem
    os.environ["ZHENG_OBSERVABLES"] = args.observables
    os.chdir(project)
    sys.path.insert(0, str(project))
    sys.path.insert(0, str(ROOT))
    import AntiGen_paths                                             # noqa: F401  (repository root on sys.path)
    import Model_run
    from Model_generate import update_antimony_model
    from Modules import Solver_settings
    from pyantigen.engine.Backend import configure_spec
    from pyantigen.engine.Model_optimize import setup_optimization_from_groups, _NO_DIAGNOSTICS

    MODEL_NAME, REPO_ROOT = Model_run.MODEL_NAME, Model_run.REPO_ROOT
    if args.no_exact:
        from pyantigen.engine.Analytic_gradient import LossGradient
        LossGradient.exact = False

    update_antimony_model()
    entry = next(iter(Model_run.EXPERIMENT_dict.values()))           # the project's one experiment
    experiment, plot = entry["EXPERIMENT"], entry["plot"]

    spec = Model_run.get_OPTIMIZATION(args.spec)
    spec.x0 = jittered_start(spec.x0, spec.bounds if isinstance(spec.bounds, dict) else None, args.seed)
    if args.max_evals:
        opts = dict((spec.optimizer_kwargs or {}).get("options") or {})
        opts.update({"maxiter": args.max_evals, **({"maxfev": args.max_evals} if spec.method == "Nelder-Mead" else {})})
        spec.optimizer_kwargs = {**(spec.optimizer_kwargs or {}), "options": opts}
    spec = configure_spec(spec, args.backend, Solver_settings, stage="fit")
    backend = "amici" if (spec.optimizer_kwargs or {}).get("analytic_gradient") else "roadrunner"

    shutil.rmtree(Path(REPO_ROOT) / "results" / MODEL_NAME / "fits", ignore_errors=True)   # never resume or reuse a fit
    settings = {"run_steady_state_first": False, "Verbose": False, "save_SBML?": False, "MODEL_NAME": MODEL_NAME,
                "REPO_ROOT": REPO_ROOT, "slice_analysis": False, "fit_mode": "optimize", "reuse_fit": False,
                "run_label": args.label}
    settings.update(_NO_DIAGNOSTICS)

    t0 = time.time()
    opt = setup_optimization_from_groups(settings, spec, {"EXPERIMENT": experiment, "plot": plot})
    wall = time.time() - t0

    rep = next(iter(experiment.replicates.values()))
    elem = next(iter(spec.groups.values()))["loss_elements"][0]
    lc = elem["loss_config"](rep) if callable(elem["loss_config"]) else elem["loss_config"]
    known = any("sigma_method" in o for o in lc["observables"])
    result = {
        "problem": args.problem, "label": args.label, "route": f"PyAntiGen / {backend} / {spec.method}", "spec": args.spec,
        "backend": backend, "method": spec.method, "sigma": "declared" if known else "profiled",
        "exact_derivatives": not args.no_exact, "observables": args.observables,
        "derivatives": ("finite differences of the whole loss" if args.no_exact else
                        "chain rule, exact" if args.observables == "column" else "chain rule, observable differenced"),
        "seed": args.seed, "note": args.note,
        "params": dict(zip(spec.param_names, [float(v) for v in np.asarray(opt["x"], dtype=float)])),
        "objective": float(opt["fun"]), "success": bool(opt["success"]), "message": str(opt["message"]),
        "nfev": opt.get("nfev"), "nit": opt.get("nit"), "wall_s": round(wall, 1),
        "platform": platform.platform(), "python": platform.python_version(),
    }
    out = HERE / args.problem / "fits" / f"{args.label}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=1), encoding="utf-8")
    print(f"wrote {out}: objective {result['objective']:.4f}, nfev {result['nfev']}, {wall:.1f} s")


if __name__ == "__main__":
    main()

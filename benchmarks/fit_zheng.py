"""Fit Zheng_PNAS2012 with one PyAntiGen route and save the result in the format compare_fits.py reads.

    python benchmarks/fit_zheng.py --label pyantigen_amici_lbfgsb --spec OPTIMIZATION_Zheng_amici

Routes (specs in problems/Zheng_PNAS2012/Modules/Optimizer_settings.py):

    OPTIMIZATION_Zheng         Nelder-Mead on RoadRunner, sigma profiled       (any platform)
    OPTIMIZATION_Zheng_amici   L-BFGS-B with the AMICI gradient, sigma profiled (WSL / Linux)
    OPTIMIZATION_Zheng_lm      Levenberg-Marquardt on AMICI, sigma declared     (WSL / Linux)

Every route starts from the same perturbed point (seed 1: each parameter times 10**U(-0.3, 0.3), PyAntiGen's parameter
order), the start of benchmarks/amici_direct_zheng.py too.  Writes benchmarks/fits/<label>.json.
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
PROJECT = ROOT / "problems" / "Zheng_PNAS2012"


def jittered_start(x0: dict, seed: int) -> dict:
    """x * 10**U(-0.3, 0.3), clipped to PEtab's bounds [1e-5, 1e3]; seed 0 keeps PEtab's nominal point."""
    if not seed:
        return dict(x0)
    rng = np.random.default_rng(seed)
    return {k: float(np.clip(v * 10 ** rng.uniform(-0.3, 0.3), 1e-5, 1e3)) for k, v in x0.items()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--label", required=True, help="name of the result file (benchmarks/fits/<label>.json)")
    ap.add_argument("--spec", default="OPTIMIZATION_Zheng", help="Optimization in Modules/Optimizer_settings.py")
    ap.add_argument("--seed", type=int, default=1, help="start perturbation (0: PEtab's nominal point)")
    ap.add_argument("--max-evals", type=int, default=None, help="evaluation budget (default: the spec's)")
    ap.add_argument("--backend", default=None, choices=["roadrunner", "amici", "auto"])
    ap.add_argument("--no-exact", action="store_true", help="finite differences of the whole loss instead of the chain rule")
    ap.add_argument("--observables", choices=["column", "callable"], default="column",
                    help="column: bare result columns; callable: the same numbers from a function of the table, which the "
                         "engine cannot take as a bare column (as a difference of simulations is not)")
    ap.add_argument("--note", default="")
    args = ap.parse_args(argv)

    os.environ["ZHENG_OBSERVABLES"] = args.observables
    os.chdir(PROJECT)
    sys.path.insert(0, str(PROJECT))
    sys.path.insert(0, str(ROOT))
    import AntiGen_paths                                             # noqa: F401  (repository root on sys.path)
    from Model_run import get_EXPERIMENT, get_OPTIMIZATION, plot_results, MODEL_NAME, REPO_ROOT
    from Model_generate import update_antimony_model
    from Modules import Solver_settings
    from pyantigen.engine.Backend import configure_spec
    from pyantigen.engine.Model_optimize import setup_optimization_from_groups, _NO_DIAGNOSTICS

    if args.no_exact:
        from pyantigen.engine.Analytic_gradient import LossGradient
        LossGradient.exact = False

    update_antimony_model()
    spec = get_OPTIMIZATION(args.spec)
    spec.x0 = jittered_start(spec.x0, args.seed)
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
    try:
        opt = setup_optimization_from_groups(settings, spec, {"EXPERIMENT": get_EXPERIMENT("EXPERIMENT_Zheng"),
                                                              "plot": plot_results})
    except RuntimeError as exc:                                       # LM that spent its budget: keep what it reached
        print("fit ended with:", exc)
        raise
    wall = time.time() - t0

    known = any("sigma_method" in o for o in spec.groups["Zheng"]["loss_elements"][0]["loss_config"]({"Label": "x"})["observables"])
    result = {
        "label": args.label, "route": f"PyAntiGen / {backend} / {spec.method}", "spec": args.spec,
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
    out = HERE / "fits" / f"{args.label}.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(result, indent=1), encoding="utf-8")
    print(f"wrote {out}: objective {result['objective']:.4f}, nfev {result['nfev']}, {wall:.1f} s")


if __name__ == "__main__":
    main()

import sys

import AntiGen_paths  # noqa: F401  (puts the repository root on sys.path before Modules is imported)

from pyantigen.engine.Model_simulate import setup_simulation
from pyantigen.engine.Model_optimize import (
    setup_optimization_from_groups,
    DIAGNOSTICS_PRESETS,
    _NO_DIAGNOSTICS,
)
from Modules.Plots import *
from Modules.Experiment import get_EXPERIMENT
from Modules.Optimizer_settings import get_OPTIMIZATION
from AntiGen_paths import MODEL_NAME, REPO_ROOT
from Model_generate import update_antimony_model

EXPERIMENT_dict = {
    "Zheng": {'EXPERIMENT': get_EXPERIMENT('EXPERIMENT_Zheng'), 'plot': plot_results, 'opt_settings_key': 'Zheng'},
}

OPTIMIZATION_REGISTRY = {
    "Zheng": {
        'opt_key': ['OPTIMIZATION_Zheng'],
        'experiment': 'EXPERIMENT_Zheng',
        'title': "Zheng_PNAS2012: fit inflowp and 44 methylation rate constants",
        'description': "All 45 parameters on a log10 scale in [1e-5, 1e3], started at PEtab's nominal values.",
    },
    "ZhengAMICI": {
        'opt_key': ['OPTIMIZATION_Zheng_amici'],
        'experiment': 'EXPERIMENT_Zheng',
        'title': "Zheng_PNAS2012: the same fit on AMICI (L-BFGS-B with analytic gradient)",
        'description': "Needs AMICI: run from WSL or Linux (python -m pyantigen.setup_env).",
    },
    "ZhengLM": {
        'opt_key': ['OPTIMIZATION_Zheng_lm'],
        'experiment': 'EXPERIMENT_Zheng',
        'title': "Zheng_PNAS2012: Levenberg-Marquardt on AMICI, sigma declared (PEtab's estimate)",
        'description': "Needs AMICI. LM is least squares, so every sigma must be known.",
    },
}

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Tellurium/RoadRunner Simulation and Optimization Runner")
    parser.add_argument("--optimize", type=str, choices=list(OPTIMIZATION_REGISTRY.keys()), help="Optimization mode")
    parser.add_argument("--simulate", type=str, choices=list(EXPERIMENT_dict.keys()), help="Simulation mode")
    parser.add_argument("--diagnostics", type=str, default="_NO_DIAGNOSTICS", choices=list(DIAGNOSTICS_PRESETS.keys()),
                        help="Diagnostics settings preset to run (default: _NO_DIAGNOSTICS)")
    parser.add_argument("--backend", type=str, default=None, choices=["auto", "roadrunner", "amici"],
                        help="Solver: roadrunner, amici (WSL/Linux; see `python -m pyantigen.setup_env`), or auto. "
                             "Overrides Solver_settings.BACKEND and PYANTIGEN_BACKEND.")
    parser.add_argument("--compare", action="store_true", dest="compare",
                        help="With --simulate --backend amici: also run RoadRunner and print the largest "
                             "difference between the two trajectories, column by column.")
    parser.add_argument("--no-fit", action="store_true", dest="no_fit",
                        help="Skip the optimizer and evaluate the spec's x0 instead, then run the requested diagnostics.")
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.error("No arguments provided.")

    if args.compare and not (args.simulate and args.backend == "amici"):
        parser.error("--compare needs --simulate with --backend amici")

    update_antimony_model()

    run_settings = {
        "run_steady_state_first": False,
        "Verbose": True,
        "save_SBML?": False,
        "MODEL_NAME": MODEL_NAME,
        "REPO_ROOT": REPO_ROOT,
    }

    if args.optimize:
        opt_info = OPTIMIZATION_REGISTRY[args.optimize]
        run_settings.update({
            "slice_analysis": False,
            "fit_mode": "evaluate_x0" if args.no_fit else "optimize",
            "run_label": args.optimize,
        })
        run_settings.update(DIAGNOSTICS_PRESETS.get(args.diagnostics, _NO_DIAGNOSTICS))
        exp_short = opt_info["experiment"].replace("EXPERIMENT_", "")
        experiment_arg = {
            "EXPERIMENT": get_EXPERIMENT(opt_info["experiment"]),
            "plot": EXPERIMENT_dict.get(exp_short, {}).get("plot"),
        }
        print(f"{args.optimize}: {opt_info['title']}\n{opt_info['description']}")
        from Modules import Solver_settings as _solver_settings
        from pyantigen.engine.Backend import BackendUnavailable, configure_spec
        for opt_key in opt_info["opt_key"]:
            try:
                spec = configure_spec(get_OPTIMIZATION(opt_key), args.backend, _solver_settings, stage="fit")
            except BackendUnavailable as exc:
                print(f"\n*** STOPPED: {exc}", flush=True)
                sys.exit(2)
            setup_optimization_from_groups(run_settings, spec, experiment_arg)

    elif args.simulate:
        # --simulate is RoadRunner unless --backend amici is given (not "auto": nothing is fitted here)
        run_settings["backend"] = args.backend if args.backend == "amici" else "roadrunner"
        run_settings["compare_backends"] = args.compare
        if args.backend == "amici":
            from Modules import Solver_settings as _solver_settings
            run_settings["amici_cache_dir"] = getattr(_solver_settings, "AMICI_CACHE_DIR", None)
            run_settings["amici_seed"] = getattr(_solver_settings, "AMICI_SEED", None)
        setup_simulation(run_settings, EXPERIMENT_dict[args.simulate])

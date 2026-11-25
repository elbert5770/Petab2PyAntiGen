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
}

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Tellurium/RoadRunner Simulation and Optimization Runner")
    parser.add_argument("--optimize", type=str, choices=list(OPTIMIZATION_REGISTRY.keys()), help="Optimization mode")
    parser.add_argument("--simulate", type=str, choices=list(EXPERIMENT_dict.keys()), help="Simulation mode")
    parser.add_argument("--diagnostics", type=str, default="_NO_DIAGNOSTICS", choices=list(DIAGNOSTICS_PRESETS.keys()),
                        help="Diagnostics settings preset to run (default: _NO_DIAGNOSTICS)")
    parser.add_argument("--no-fit", action="store_true", dest="no_fit",
                        help="Skip the optimizer and evaluate the spec's x0 instead, then run the requested diagnostics.")
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.error("No arguments provided.")

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
        for opt_key in opt_info["opt_key"]:
            setup_optimization_from_groups(run_settings, get_OPTIMIZATION(opt_key), experiment_arg)

    elif args.simulate:
        setup_simulation(run_settings, EXPERIMENT_dict[args.simulate])

import os

import pandas as pd

from AntiGen_paths import REPO_ROOT

# the last measurement (the receptor starts in state R at time 0, PEtab's time zero)
T_END = float(pd.read_csv(os.path.join(REPO_ROOT, "data", "Sneyd_PNAS2002.csv"))["time"].max())

# Solver backend: "roadrunner", "amici" or "auto" (the stage default).  Overridden by
# `Model_run.py --backend` and by the PYANTIGEN_BACKEND environment variable.
BACKEND = "auto"
# Where AMICI keeps compiled models.  None: PYANTIGEN_AMICI_CACHE, else ~/.cache/pyantigen/amici_models.
AMICI_CACHE_DIR = None
# How a data-floored block gets its sigma ("auto" picks "fixed" for AMICI multistart/profile/FIM).
# Sneyd declares no floors, so this changes nothing here.
SIGMA_FLOOR = "auto"


def make_solver_settings(blocks, abs_tol=1e-10, rel_tol=1e-10, max_steps=200000):
    """Factory for solver settings dicts."""
    return {
        'integrator': 'cvode',
        'absolute_tolerance': abs_tol,
        'relative_tolerance': rel_tol,
        'maximum_num_steps': max_steps,
        'stiff': True,
        'variable_step_size': True,
        'simulation_blocks': blocks,
    }


def solver_settings_Sneyd(replicate):
    # one output row every 0.01 s: the measurement times (multiples of 0.07) are among them
    return make_solver_settings({'observe': {'start': 0, 'end': T_END, 'n_points': round(T_END / 0.01) + 1}})

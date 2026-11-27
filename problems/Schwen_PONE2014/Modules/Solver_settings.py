import math
import os

import pandas as pd

from AntiGen_paths import REPO_ROOT

# Solver backend: "roadrunner", "amici" or "auto" (the stage default).  Overridden by
# `Model_run.py --backend` and by the PYANTIGEN_BACKEND environment variable.
BACKEND = "auto"
# Where AMICI keeps compiled models.  None: PYANTIGEN_AMICI_CACHE, else ~/.cache/pyantigen/amici_models.
AMICI_CACHE_DIR = None
# How a data-floored block gets its sigma ("auto" picks "fixed" for AMICI multistart/profile/FIM). Schwen declares no floors,
# so this changes nothing here.
SIGMA_FLOOR = "auto"

_DATA = os.path.join(REPO_ROOT, "data", "Schwen_PONE2014.csv")
_UNIT = 0.05                      # every measurement time is a multiple of this (minutes)


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


def solver_settings_Schwen(replicate):
    """Simulate from 0 to this condition's last measurement, on a grid that holds every measurement time (about 100 points).
    A condition measured only at t = 0 (the first of the cytometry doses) is simulated to one minute."""
    df = pd.read_csv(_DATA)
    times = df.loc[df["condition"] == replicate["Label"], "time"].to_numpy()
    units = sorted({round(t / _UNIT) for t in times} | {0})
    if units[-1] == 0:
        units.append(round(1 / _UNIT))
    t_end = units[-1]
    g = math.gcd(*[u for u in units if u > 0])                  # the coarsest grid that holds all of them
    steps = math.ceil(100 / (t_end // g)) * (t_end // g)         # subdivide to get at least 100 intervals
    return make_solver_settings({'observe': {'start': 0, 'end': t_end * _UNIT, 'n_points': steps + 1}})

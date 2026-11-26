# The experiment: the cell is pre-equilibrated with every input off (PEtab's model1_data1), then one input is switched on at
# once and a second 24 h later (Events.py), and the signals are followed for 30 h.  Two simulation blocks and the events
# stand for PEtab's "preequilibration -> simulation condition".

# Pre-equilibration: the cascade settles within ~1e4 h at every parameter set tried (the pre-equilibrated state agrees with a
# steady-state solve to 1e-13); 1e5 h is a 10x margin.  Untracked: only the state it leaves matters.
T_PREEQUIL = 1e5
T_OBSERVE = 30.0          # last measurement, h after the switch
N_OBSERVE = 3001          # 0.01 h: contains every measurement time (24.08, 24.17, ...)


# Solver backend: "roadrunner", "amici" or "auto" (the stage default).  Overridden by
# `Model_run.py --backend` and by the PYANTIGEN_BACKEND environment variable.
BACKEND = "auto"
# Where AMICI keeps compiled models.  None: PYANTIGEN_AMICI_CACHE, else ~/.cache/pyantigen/amici_models.
AMICI_CACHE_DIR = None
# How a data-floored block gets its sigma ("auto" picks "fixed" for AMICI multistart/profile/FIM).
# Weber declares no floors, so this changes nothing here.
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


def solver_settings_Weber(replicate):
    return make_solver_settings({
        'preequil': {'start': 0, 'end': T_PREEQUIL, 'n_points': 2, 'tracked': False},
        'observe': {'start': T_PREEQUIL, 'end': T_PREEQUIL + T_OBSERVE, 'n_points': N_OBSERVE},
    })

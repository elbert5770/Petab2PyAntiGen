# The experiment: the culture is pre-equilibrated with dilution on (dilution = 1,
# inflow = inflowp), then dilution is switched off and the 15 H3 states are followed
# for 25 h.  PEtab says "preequilibration -> simCondition"; here the same thing is
# two simulation blocks and one event (Events.py).

# Pre-equilibration: every state is diluted at rate inflowp (~0.03 /h at the fit), so
# 1e5 h is > 3000 relaxation times.  Untracked: only the state it leaves matters.
T_PREEQUIL = 1e5
T_OBSERVE = 25.0          # last measurement, h after the switch


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


def solver_settings_Zheng(replicate):
    return make_solver_settings({
        'preequil': {'start': 0, 'end': T_PREEQUIL, 'n_points': 2, 'tracked': False},
        'observe': {'start': T_PREEQUIL, 'end': T_PREEQUIL + T_OBSERVE, 'n_points': 501},
    })

"""Experiment registry: one 'replicate' entry per PEtab simulation condition."""
from .Data import *
from .Events import *
from .Loss_config import *
from .Observed_species import *
from .Solver_settings import *
from .Update_parameters import *
from .Update_opt_parameters import *

from dataclasses import dataclass, field


@dataclass
class Experiment:
    # Which replicates are fitted together, and how their losses are weighted, is a property of the optimization:
    # Optimization.groups in Optimizer_settings.py. The experiment only says what can be simulated.
    replicates: dict = field(default_factory=dict)


def make_replicate(config, **meta):
    """Build one replicate entry from a config dict (see PyAntiGen's template Experiment.py)."""
    return {**dict(config), **meta}


def _build_experiment():
    exp = Experiment()
    for row in load_conditions().to_dict("records"):          # 19 PEtab conditions
        key = row["condition"]
        exp.replicates[key] = make_replicate({
            "Label":                 key,
            "Events":                generate_no_events,
            "Data":                  load_Schwen_data,
            "Observed_species":      all_species,
            "Solver_settings":       solver_settings_Schwen,
            "Update_parameters":     update_Schwen,
            "Update_opt_parameters": update_opt_Schwen,
            "no_inputs":             True,      # nothing is clock-driven: AMICI needs no input schedule
        }, **row)
    return exp


EXPERIMENT_Schwen = _build_experiment()


def get_EXPERIMENTS():
    """Returns a dict mapping experiment name -> Experiment."""
    return {k: v for k, v in globals().items() if k.startswith('EXPERIMENT_')}


def get_EXPERIMENT(name):
    """Returns a single Experiment by name."""
    return globals().get(name)

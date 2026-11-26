from dataclasses import dataclass, field

from .Data import load_conditions
from .Loss_config import *


@dataclass
class Optimization:
    param_names: list
    x0: list
    bounds: list = None
    method: str = "Nelder-Mead"
    optimizer_kwargs: dict = field(default_factory=dict)
    group_normalization: str = "mean_over_groups"  # "mean_over_groups" | "sum_over_groups"
    groups: dict = field(default_factory=dict)       # nested group/loss configuration
    passive_simulations: list = field(default_factory=list) # passive simulations to run for plotting
    # None/"lin" keeps a parameter linear; "log10" fits log10(p). A string applies to all.
    # x0 and bounds stay in linear units; so does everything reported back.
    parameter_scale: object = None


# PEtab's bounds: R in [1e-3, 1e2]; the two rates in [1e-8, 1]
BOUNDS = {"R": (1e-3, 1e2), "gamma_1": (1e-8, 1.0), "gamma_2": (1e-8, 1.0)}


def _build_Zhao_fit(method="Nelder-Mead", optimizer_kwargs=None):
    """Fit R, gamma_1 and gamma_2 of every (region, stage): 21 parameters, each condition's three on their own data."""
    conditions = load_conditions().to_dict("records")
    names, x0, bounds = [], {}, {}
    for row in conditions:
        for p in BOUNDS:
            name = f"{p}_{row['stage']}_{row['region']}"          # R_Stage_I_Wuhan
            names.append(name)
            x0[name], bounds[name] = row[p], BOUNDS[p]
    return Optimization(
        param_names=names, x0=x0, bounds=bounds,
        method=method,
        optimizer_kwargs=optimizer_kwargs or {"options": {"maxiter": 20000, "maxfev": 20000}},
        group_normalization="sum_over_groups",
        parameter_scale="log10",
        groups={"Zhao": {"group_weight": 1.0, "loss_elements": [
            {"simulation": f"{row['region']}_{row['stage']}", "loss_config": Zhao_loss_config, "weight": 1.0}
            for row in conditions]}},
        passive_simulations=[],
    )


OPTIMIZATION_Zhao = _build_Zhao_fit()

# The same fit on AMICI: analytic_gradient makes AMICI the backend and L-BFGS-B uses its gradient. Each fitted name reaches the
# model through Update_opt_parameters; AMICI's setup replays that hook to learn which model parameter it moves.
OPTIMIZATION_Zhao_amici = _build_Zhao_fit(
    method="L-BFGS-B", optimizer_kwargs={"analytic_gradient": True, "options": {"maxiter": 2000}})


def get_OPTIMIZATION(name):
    """Returns a single Optimization configuration by name."""
    return globals().get(name)

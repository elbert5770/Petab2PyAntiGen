import os
from dataclasses import dataclass, field

import pandas as pd

from antimony_modules.Calcium.ip3r_Sneyd_PNAS2002 import RATE_CONSTANTS
from .Data import load_conditions
from .Loss_config import *
from AntiGen_paths import MODEL_NAME, REPO_ROOT


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


def _nominal():
    """Start values: the model's own parameter values (antimony_models/<Model>/..._parameters.csv)."""
    df = pd.read_csv(os.path.join(REPO_ROOT, "antimony_models", MODEL_NAME, f"{MODEL_NAME}_parameters.csv"))
    return dict(zip(df["Parameter"], df["Value"].astype(float)))


def _build_Sneyd_fit(method="Nelder-Mead", optimizer_kwargs=None, loss_config=None):
    """Fit the 14 rate constants, all on a log10 scale in [1e-3, 1e5] (PEtab's parameters.tsv); sigma is profiled out."""
    nominal = _nominal()
    return Optimization(
        param_names=list(RATE_CONSTANTS),
        x0={n: nominal[n] for n in RATE_CONSTANTS},
        bounds={n: (1e-3, 1e5) for n in RATE_CONSTANTS},
        method=method,
        optimizer_kwargs=optimizer_kwargs or {"options": {"maxiter": 20000, "maxfev": 20000}},
        group_normalization="sum_over_groups",
        parameter_scale="log10",
        groups={"Sneyd": {"group_weight": 1.0, "loss_elements": [
            {"simulation": row["condition"], "loss_config": loss_config or Sneyd_loss_config, "weight": 1.0,
             "sigma_block": "Sneyd"}                  # one sigma for all nine conditions
            for row in load_conditions().to_dict("records")]}},
        passive_simulations=[],
    )


OPTIMIZATION_Sneyd = _build_Sneyd_fit()

# The same fit on AMICI: analytic_gradient makes AMICI the backend, and L-BFGS-B uses its gradient.
OPTIMIZATION_Sneyd_amici = _build_Sneyd_fit(
    method="L-BFGS-B", optimizer_kwargs={"analytic_gradient": True, "options": {"maxiter": 2000}})

# Levenberg-Marquardt on AMICI: least squares needs a known sigma, so this one declares PEtab's.
OPTIMIZATION_Sneyd_lm = _build_Sneyd_fit(
    method="LM", optimizer_kwargs={"analytic_gradient": True, "options": {"maxiter": 1000}},
    loss_config=Sneyd_loss_config_known_sigma)


def get_OPTIMIZATION(name):
    """Returns a single Optimization configuration by name."""
    return globals().get(name)

import os
from dataclasses import dataclass, field

import pandas as pd

from antimony_modules.Histone.methylation_Zheng_PNAS2012 import rate_constants
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


def _build_Zheng_fit(method="Nelder-Mead", optimizer_kwargs=None, loss_config=None):
    """Fit inflowp and the 44 rate constants, all on a log10 scale in [1e-5, 1e3] (PEtab's parameters.tsv)."""
    names = ["inflowp"] + rate_constants()
    nominal = _nominal()
    return Optimization(
        param_names=names,
        x0={n: nominal[n] for n in names},
        bounds={n: (1e-5, 1e3) for n in names},
        method=method,
        optimizer_kwargs=optimizer_kwargs or {"options": {"maxiter": 20000, "maxfev": 20000}},
        group_normalization="mean_over_groups",
        parameter_scale="log10",
        groups={
            "Zheng": {
                "group_weight": 1.0,
                "loss_elements": [
                    {"simulation": "Zheng_PNAS2012", "loss_config": loss_config or Zheng_loss_config, "weight": 1.0},
                ],
            }
        },
        passive_simulations=[],
    )


OPTIMIZATION_Zheng = _build_Zheng_fit()

# The same fit on AMICI: analytic_gradient makes AMICI the backend, and L-BFGS-B uses its gradient.
# About 400 evaluations reach the optimum (Nelder-Mead on RoadRunner needs ~20,000).
OPTIMIZATION_Zheng_amici = _build_Zheng_fit(
    method="L-BFGS-B",
    optimizer_kwargs={"analytic_gradient": True, "options": {"maxiter": 1000}})


# Levenberg-Marquardt on AMICI: least squares needs known sigmas, so this one declares PEtab's estimate of sigma.
OPTIMIZATION_Zheng_lm = _build_Zheng_fit(
    method="LM", optimizer_kwargs={"analytic_gradient": True, "options": {"maxiter": 1000}},
    loss_config=Zheng_loss_config_known_sigma)


def get_OPTIMIZATION(name):
    """Returns a single Optimization configuration by name."""
    return globals().get(name)

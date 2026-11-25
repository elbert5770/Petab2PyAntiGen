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


def _build_Zheng_fit():
    """Fit inflowp and the 44 rate constants, all on a log10 scale in [1e-5, 1e3] (PEtab's parameters.tsv)."""
    names = ["inflowp"] + rate_constants()
    nominal = _nominal()
    return Optimization(
        param_names=names,
        x0={n: nominal[n] for n in names},
        bounds={n: (1e-5, 1e3) for n in names},
        method="Nelder-Mead",
        optimizer_kwargs={"options": {"maxiter": 20000, "maxfev": 20000}},
        group_normalization="mean_over_groups",
        parameter_scale="log10",
        groups={
            "Zheng": {
                "group_weight": 1.0,
                "loss_elements": [
                    {"simulation": "Zheng_PNAS2012", "loss_config": Zheng_loss_config, "weight": 1.0},
                ],
            }
        },
        passive_simulations=[],
    )


OPTIMIZATION_Zheng = _build_Zheng_fit()


def get_OPTIMIZATION(name):
    """Returns a single Optimization configuration by name."""
    return globals().get(name)

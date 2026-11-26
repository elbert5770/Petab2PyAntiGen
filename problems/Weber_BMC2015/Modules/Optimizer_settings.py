import os
from dataclasses import dataclass, field

import pandas as pd

from .Events import CONDITIONS
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


# PEtab's bounds, by parameter family (parameters.tsv); every one is fitted on a log10 scale
BOUNDS = {"a": (1e-4, 1e2),             # first-order degradation
          "m": (1e2, 1e10),             # saturation constants
          "p": (1e-5, 1e5),             # activation / transfer rates
          "pu": (1e-1, 1e8),            # input strengths (pu2, a rate that is always off, is not estimated)
          "s": (1e2, 1e7),              # constant synthesis
          "scale": (1e-5, 1e3)}         # observable scale factors


def _nominal():
    """Start values: the model's own parameter values (antimony_models/<Model>/..._parameters.csv)."""
    df = pd.read_csv(os.path.join(REPO_ROOT, "antimony_models", MODEL_NAME, f"{MODEL_NAME}_parameters.csv"))
    return dict(zip(df["Parameter"], df["Value"].astype(float)))


def _family(name):
    return name.split("_")[0].rstrip("0123456789")


def _build_Weber_fit(method="Nelder-Mead", optimizer_kwargs=None):
    """Fit the 31 parameters PEtab estimates besides its five noise parameters, which are profiled out."""
    nominal = _nominal()
    names = [n for n in nominal if _family(n) in BOUNDS]
    return Optimization(
        param_names=names,
        x0={n: nominal[n] for n in names},
        bounds={n: BOUNDS[_family(n)] for n in names},
        method=method,
        optimizer_kwargs=optimizer_kwargs or {"options": {"maxiter": 20000, "maxfev": 20000}},
        group_normalization="sum_over_groups",
        parameter_scale="log10",
        groups={"Weber": {"group_weight": 1.0, "loss_elements": [
            {"simulation": key, "loss_config": Weber_loss_config, "weight": 1.0, "sigma_block": "Weber"} for key in CONDITIONS]}},
        passive_simulations=[],
    )


OPTIMIZATION_Weber = _build_Weber_fit()

# The same fit on AMICI: analytic_gradient makes AMICI the backend, and L-BFGS-B uses its gradient.
OPTIMIZATION_Weber_amici = _build_Weber_fit(
    method="L-BFGS-B", optimizer_kwargs={"analytic_gradient": {"rtol": 1e-8, "atol": 1e-10, "retry": [[1e-6, 1e-8]]},
                                         "options": {"maxiter": 2000}})


def get_OPTIMIZATION(name):
    """Returns a single Optimization configuration by name."""
    return globals().get(name)

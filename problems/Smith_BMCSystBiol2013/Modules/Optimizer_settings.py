import os
from dataclasses import dataclass, field

import pandas as pd

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


# What PEtab estimates (parameters.tsv): 16 rate constants, bounded by 1/100 and 100 times the nominal value, and the 9
# observable scales, bounded by 1e-4 and 1e4. Every other parameter is fixed. All on a log10 scale; there is no noise parameter.
RATE_CONSTANTS = ("k1", "kminus1", "k30f", "k30r", "k35f", "k7", "kminus7a", "kminus7b", "k8", "kminus8",
                  "k42f", "k42r", "k43f", "k43r", "k32f", "k32r")
SCALES = ("sc_PI3K", "sc_GLUT_2B", "sc_BINS", "sc_PIRS", "sc_PTP", "sc_GLUT_3B_120", "sc_GLUT_3B_240", "sc_SOD2", "sc_FOXO1")


def _nominal():
    """Start values: the model's own parameter values (antimony_models/<Model>/..._parameters.csv)."""
    df = pd.read_csv(os.path.join(REPO_ROOT, "antimony_models", MODEL_NAME, f"{MODEL_NAME}_parameters.csv"))
    return dict(zip(df["Parameter"], df["Value"].astype(float)))


def _build_Smith_fit(method="Nelder-Mead", optimizer_kwargs=None):
    """Fit PEtab's 25 parameters to the 62 measurements of 36 conditions, all sigmas declared (1)."""
    nominal = _nominal()
    names = list(RATE_CONSTANTS + SCALES)
    bounds = {n: (nominal[n] / 100, nominal[n] * 100) if n in RATE_CONSTANTS else (1e-4, 1e4) for n in names}
    return Optimization(
        param_names=names,
        x0={n: nominal[n] for n in names},
        bounds=bounds,
        method=method,
        optimizer_kwargs=optimizer_kwargs or {"options": {"maxiter": 20000, "maxfev": 20000}},
        group_normalization="sum_over_groups",
        parameter_scale="log10",
        groups={"Smith": {"group_weight": 1.0, "loss_elements": [
            {"simulation": row["condition"], "loss_config": Smith_loss_config, "weight": 1.0}
            for row in load_conditions().to_dict("records")]}},
        passive_simulations=[],
    )


OPTIMIZATION_Smith = _build_Smith_fit()

# The same fit on AMICI: analytic_gradient makes AMICI the backend, and L-BFGS-B uses its gradient.
# Akt and PKC switch on through piecewise(PIP3 - threshold, PIP3 > threshold, 0), which AMICI integrates as an event. Where PIP3
# crosses the threshold CVODE can return "a root after reinitialization" (status -99), whether it does is knife-edge in the
# parameters (a 1e-16 change decides it), and it is worse at tight tolerances: a looser tolerance, tried in turn, succeeds.
GRADIENT = {"rtol": 1e-8, "atol": 1e-10, "retry": [[1e-6, 1e-8], [1e-5, 1e-7], [1e-4, 1e-6]]}
OPTIMIZATION_Smith_amici = _build_Smith_fit(
    method="L-BFGS-B", optimizer_kwargs={"analytic_gradient": GRADIENT, "options": {"maxiter": 2000}})

# Levenberg-Marquardt on AMICI (every sigma is declared, as least squares needs). The default projected-gradient tolerance is 0.01 nats;
# at an objective of ~2.5e4 whose integration noise is ~1e-4 relative, steps from a gradient below ~2 are rejected (the predicted
# decrease is smaller than the noise), so a fit that has stopped improving is reported as not converged and discarded. 2 nats is
# 1e-4 of the objective.
OPTIMIZATION_Smith_lm = _build_Smith_fit(
    method="LM", optimizer_kwargs={"analytic_gradient": GRADIENT, "options": {"maxiter": 1000, "pg_tol": 2.0}})


def get_OPTIMIZATION(name):
    """Returns a single Optimization configuration by name."""
    return globals().get(name)

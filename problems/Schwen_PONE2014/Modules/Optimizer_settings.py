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
    # None/"lin" keeps a parameter linear; "log10" fits log10(p). A string applies to all, a dict names them.
    # x0 and bounds stay in linear units; so does everything reported back.
    parameter_scale: object = None


# What PEtab estimates (parameters.tsv), besides the two noise parameters that are profiled out: 28 parameters. Bounds by family.
RATES = ("ka1", "ka2fold", "kd1", "kd2fold", "kin", "kin2", "kon_unspec", "koff_unspec", "kout", "kout2", "kout_frag",
         "ini_R1", "ini_R2fold", "offset", "scale")
BOUNDS = {**{p: (1e-5, 1e3) for p in RATES},
          **{f"offset_nExpID{k}": (1e-5, 1e5) for k in (1, 2, 3, 4)},
          **{f"km_nExpID{k}": (1e3, 1e8) for k in (1, 2, 3, 4)},
          **{f"scaleElisa_nExpID{k}": (0.1, 1.0) for k in (1, 2, 3, 4)},          # linear in PEtab
          "fragments": (0.0, 1.0)}                                                 # linear in PEtab
LINEAR = ("fragments", "scaleElisa_nExpID1", "scaleElisa_nExpID2", "scaleElisa_nExpID3", "scaleElisa_nExpID4")

# PEtab's priors (parameterScaleNormal, on the log10 scale): the mean in decades and sd = 3 decades, on the six rate parameters that
# have one. On AMICI they are log-normal priors of the engine (``analytic_gradient["priors"]``); RoadRunner's objective has none.
PRIOR_MEAN = {"ini_R1": 1.809765366817905, "ini_R2fold": 1.196922192866968, "ka1": -4.006774320046476,
              "ka2fold": -0.018613892336092, "kd1": -0.917792647716124, "kd2fold": -0.105608820817629}
PRIORS = [{"param": p, "ref": 10 ** m, "sigma": 3.0, "side": "both"} for p, m in PRIOR_MEAN.items()]


def _nominal():
    """Start values: the model's own parameter values (antimony_models/<Model>/..._parameters.csv)."""
    df = pd.read_csv(os.path.join(REPO_ROOT, "antimony_models", MODEL_NAME, f"{MODEL_NAME}_parameters.csv"))
    return dict(zip(df["Parameter"], df["Value"].astype(float)))


def _build_Schwen_fit(method="Nelder-Mead", optimizer_kwargs=None):
    """Fit PEtab's 28 parameters (and two profiled sigmas) to the 286 measurements of 19 conditions."""
    nominal = _nominal()
    names = list(BOUNDS)
    return Optimization(
        param_names=names,
        x0={n: nominal[n] for n in names},
        bounds=dict(BOUNDS),
        method=method,
        optimizer_kwargs=optimizer_kwargs or {"options": {"maxiter": 20000, "maxfev": 20000}},
        group_normalization="sum_over_groups",
        parameter_scale={n: "lin" if n in LINEAR else "log10" for n in names},
        groups={"Schwen": {"group_weight": 1.0, "loss_elements": [
            {"simulation": row["condition"], "loss_config": Schwen_loss_config, "weight": 1.0,
             "sigma_block": block_of(row["condition"])}
            for row in load_conditions().to_dict("records")]}},
        passive_simulations=[],
    )


OPTIMIZATION_Schwen = _build_Schwen_fit()

# The same fit on AMICI: analytic_gradient makes AMICI the backend, and L-BFGS-B uses its gradient. The priors are part of it.
GRADIENT = {"rtol": 1e-8, "atol": 1e-10, "retry": [[1e-6, 1e-8], [1e-5, 1e-7], [1e-4, 1e-6]], "priors": PRIORS}
OPTIMIZATION_Schwen_amici = _build_Schwen_fit(
    method="L-BFGS-B", optimizer_kwargs={"analytic_gradient": GRADIENT, "options": {"maxiter": 2000}})


def get_OPTIMIZATION(name):
    """Returns a single Optimization configuration by name."""
    return globals().get(name)

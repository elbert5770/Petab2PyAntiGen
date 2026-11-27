import os

import pandas as pd

from AntiGen_paths import REPO_ROOT

# PEtab's noise is the constant 1.0 for every observable (observables.tsv: noiseFormula = 1.0), and no noise parameter is
# estimated: the sigmas are declared, so the loss is half the sum of squared residuals plus a constant.
SIGMA = 1.0

_DATA = os.path.join(REPO_ROOT, "data", "Smith_BMCSystBiol2013.csv")


def measured(label):
    """The observables a condition has data for: the columns of its rows that are not blank."""
    df = pd.read_csv(_DATA)
    df = df[df["condition"] == label].drop(columns=["condition", "time"])
    return [c for c in df.columns if df[c].notna().any()]


def no_optimization(replicate):
    return {}


def Smith_loss_config(replicate):
    """PEtab's observables, each a rule of the model (a fitted scale times a quantity), normal noise with sigma = 1."""
    label = replicate["Label"]
    return {"observables": [{
        "observed_variable": name,
        "data_column":       name,
        "time_column":       "time",
        "data_dict_key":     label,
        "sigma_method":      "fixed",
        "sigma_value":       SIGMA,
    } for name in measured(label)]}

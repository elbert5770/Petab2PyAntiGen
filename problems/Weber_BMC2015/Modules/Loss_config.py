import os

import pandas as pd

from AntiGen_paths import REPO_ROOT

# PEtab does not estimate the noise of the three "total" observables: parameters.tsv fixes std_yPKDt, std_yPI4K3Bt and std_yCERTt
# (estimate = 0). They are declared sigmas, not profiled ones. The other five observables each estimate one sigma
# (std_<observable>), shared by every condition that measures the observable: both loss elements carry the same "sigma_block"
# (Optimizer_settings.py), so the block of an observable is one across conditions (yPKDpN0 is measured in both).
KNOWN_SIGMA = {"yPKDt": 1910527.835, "yPI4K3Bt": 887676.131, "yCERTt": 244727220.419}

_DATA = os.path.join(REPO_ROOT, "data", "Weber_BMC2015.csv")


def measured(label):
    """The observables a condition has data for: the columns of its rows that are not blank."""
    df = pd.read_csv(_DATA)
    df = df[df["condition"] == label].drop(columns=["condition", "time"])
    return [c for c in df.columns if df[c].notna().any()]


def no_optimization(replicate):
    return {}


def Weber_loss_config(replicate):
    """PEtab's eight observables: each is a result column (an assignment rule of the model), normal noise."""
    label = replicate["Label"]
    return {"observables": [{
        "observed_variable": name,
        "data_column":       name,
        "time_column":       "time",
        "data_dict_key":     label,
        **({"sigma_method": "fixed", "sigma_value": KNOWN_SIGMA[name]} if name in KNOWN_SIGMA else {}),
    } for name in measured(label)]}

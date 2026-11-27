import functools
import os

import pandas as pd

from AntiGen_paths import REPO_ROOT

# PEtab estimates two noise parameters, both on the log10 scale of the observables: IR_obs_std is shared by the three cytometry
# observables (IR1, IR2, IRsum) in every condition, std by the four ELISA series. Each is one profiled sigma: the observables
# of a group carry the same "sigma_pool" (they differ inside a condition) and its loss elements the same "sigma_block" (they differ
# across conditions); Optimizer_settings.py sets the block, here the pool.
POOL = {"observable_IR1": "IR_obs_std", "observable_IR2": "IR_obs_std", "observable_IRsum": "IR_obs_std"}   # the rest: "std"

_DATA = os.path.join(REPO_ROOT, "data", "Schwen_PONE2014.csv")


@functools.lru_cache(maxsize=None)          # the engine asks for it twice per simulation per gradient: read the file once
def measured(label):
    """The observables a condition has data for: the columns of its rows that are not blank."""
    df = pd.read_csv(_DATA)
    df = df[df["condition"] == label].drop(columns=["condition", "time"])
    return tuple(c for c in df.columns if df[c].notna().any())


def block_of(label):
    """Which noise parameter a condition's measurements share: IR_obs_std (cytometry) or std (ELISA)."""
    return POOL.get(measured(label)[0], "std")


def no_optimization(replicate):
    return {}


def Schwen_loss_config(replicate):
    """PEtab's observables, each a log10 rule of the model with normal noise on that scale."""
    label = replicate["Label"]
    return {"observables": [{
        "observed_variable": name,
        "data_column":       name,
        "time_column":       "time",
        "data_dict_key":     label,
        "sigma_pool":        POOL.get(name, "std"),
    } for name in measured(label)]}

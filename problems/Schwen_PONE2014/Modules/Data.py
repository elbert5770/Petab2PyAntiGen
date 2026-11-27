import os

import pandas as pd

from AntiGen_paths import REPO_ROOT

DATA = os.path.join(REPO_ROOT, "data")


def load_conditions():
    """One row per PEtab simulation condition: the insulin (nM) the experiment starts with."""
    return pd.read_csv(os.path.join(DATA, "Schwen_PONE2014_conditions.csv"))


def load_Schwen_data(replicate, data_path):
    """This condition's measurements, as log10 (PEtab's observables are log10-transformed): one row per time, one column per
    observable it measures (the others are blank)."""
    df = pd.read_csv(os.path.join(data_path, "Schwen_PONE2014.csv"))
    df = df[df["condition"] == replicate["Label"]].drop(columns="condition").dropna(axis=1, how="all")
    return {replicate["Label"]: df.reset_index(drop=True)}

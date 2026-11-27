import os

import pandas as pd

from AntiGen_paths import REPO_ROOT

DATA = os.path.join(REPO_ROOT, "data")


def load_conditions():
    """One row per PEtab simulation condition: the insulin dose, the H2O2 outside the cell, the two on/off switches and E2F1."""
    return pd.read_csv(os.path.join(DATA, "Smith_BMCSystBiol2013_conditions.csv"))


def load_Smith_data(replicate, data_path):
    """This condition's measurements: one row per time, one column per observable it measures (the others are blank)."""
    df = pd.read_csv(os.path.join(data_path, "Smith_BMCSystBiol2013.csv"))
    df = df[df["condition"] == replicate["Label"]].drop(columns="condition").dropna(axis=1, how="all")
    return {replicate["Label"]: df.reset_index(drop=True)}

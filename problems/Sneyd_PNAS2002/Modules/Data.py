import os

import pandas as pd

from AntiGen_paths import REPO_ROOT

DATA = os.path.join(REPO_ROOT, "data")


def load_conditions():
    """One row per condition: the constant IP3 and Ca2+ (uM) of a dose-response experiment."""
    return pd.read_csv(os.path.join(DATA, "Sneyd_PNAS2002_conditions.csv"))


def load_Sneyd_data(replicate, data_path):
    """Open probability against time, of this replicate's condition."""
    df = pd.read_csv(os.path.join(data_path, "Sneyd_PNAS2002.csv"))
    return {replicate["Label"]: df[df["condition"] == replicate["Label"]].reset_index(drop=True)}

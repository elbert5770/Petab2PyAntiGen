import os

import pandas as pd

from AntiGen_paths import REPO_ROOT

DATA = os.path.join(REPO_ROOT, "data")


def load_conditions():
    """One row per (region, stage): population, initial infections and PEtab's nominal R, gamma_1, gamma_2."""
    return pd.read_csv(os.path.join(DATA, "Zhao_QuantBiol2020_conditions.csv"))


def load_no_data(replicate, data_path):
    return {}


def load_Zhao_data(replicate, data_path):
    """Cumulative confirmed infected, by day, of this replicate's (region, stage)."""
    df = pd.read_csv(os.path.join(data_path, "Zhao_QuantBiol2020.csv"))
    return {replicate["Label"]: df[df["condition"] == replicate["condition"]].reset_index(drop=True)}

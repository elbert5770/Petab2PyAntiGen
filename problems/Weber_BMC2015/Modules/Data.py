import os

import pandas as pd

from .Solver_settings import T_PREEQUIL


def load_Weber_data(replicate, data_path):
    """This condition's measurements: one row per (time, replicate), one column per observable, blank = not measured.

    Simulation time starts at the pre-equilibration, so the measurement times are shifted by T_PREEQUIL.
    """
    df = pd.read_csv(os.path.join(data_path, 'Weber_BMC2015.csv'))
    df = df[df['condition'] == replicate["Label"]].drop(columns='condition').dropna(axis=1, how='all')
    df['time'] = df['time'] + T_PREEQUIL
    return {replicate["Label"]: df.reset_index(drop=True)}

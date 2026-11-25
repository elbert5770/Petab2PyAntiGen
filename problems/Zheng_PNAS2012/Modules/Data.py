import os

import pandas as pd

from .Solver_settings import T_PREEQUIL


def load_no_data(replicate, data_path):
    return {}


def load_Zheng_data(replicate, data_path):
    """Fractions of total H3 for the 15 states at t = 0, 5, 10, 25 h after the switch.

    One row per time, one column per observable (``observable_K27me<a>K36me<b>``).
    Simulation time starts at the pre-equilibration, so the measurement times are
    shifted by T_PREEQUIL.
    """
    df = pd.read_csv(os.path.join(data_path, 'Zheng_PNAS2012.csv'))
    df['time'] = df['time'] + T_PREEQUIL
    return {replicate["Label"]: df}

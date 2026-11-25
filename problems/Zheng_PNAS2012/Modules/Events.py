from .Solver_settings import T_PREEQUIL


def generate_no_events(replicate, df_dict):
    return ''


def Zheng_event(replicate, df_dict):
    # simCondition: dilution 1 -> 0.  Dilution is the only thing that differs from
    # the pre-equilibration condition, so the whole conditions.tsv is this one line.
    return f'at (time >= {T_PREEQUIL}): dilution = 0'

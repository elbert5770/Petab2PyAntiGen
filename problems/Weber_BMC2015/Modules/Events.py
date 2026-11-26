from .Solver_settings import T_PREEQUIL

# What each PEtab condition does after the pre-equilibration, which has every input off: (hours after the switch, input, level).
# PEtab says this with six condition parameters per condition (a flag, a dose and a time for each of two stimuli) that SBML
# rules turn into the inputs u3..u6; here it is the whole of conditions.tsv.
CONDITIONS = {
    "model1_data2": [(0, "u3", 1), (24, "u5", 1)],      # ectopic PI4K3beta on at once; PdBu (dose 1) at 24 h
    "model1_data3": [(0, "u4", 1), (24, "u6", 10)],     # ectopic CERT on at once; PKD inhibitor NB142-70 (dose 10) at 24 h
}


def generate_no_events(replicate, df_dict):
    return ''


def Weber_events(replicate, df_dict):
    return "\n".join(f"at (time >= {T_PREEQUIL + hours!r}): {name} = {level}"
                     for hours, name, level in CONDITIONS[replicate["Label"]])

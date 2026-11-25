from antimony_modules.Histone.methylation_Zheng_PNAS2012 import STATES, state

# True: one sigma for all 15 observables (PEtab's Zheng problem: a single estimated
# parameter `sigma` in every noiseFormula).  False: one sigma per observable.
SHARED_SIGMA = True


def no_optimization(replicate):
    return {}


def Zheng_loss_config(replicate):
    """The 15 observables of PEtab's observables.tsv: each state as a fraction of total H3, normal noise."""
    label = replicate["Label"]
    return {"observables": [{
        "observed_variable": f"observable_{state(a, b)}",
        "data_column":       f"observable_{state(a, b)}",
        "time_column":       "time",
        "data_dict_key":     label,
        **({"sigma_pool": "sigma"} if SHARED_SIGMA else {}),
    } for a, b in STATES]}

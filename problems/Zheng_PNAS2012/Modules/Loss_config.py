from antimony_modules.Histone.methylation_Zheng_PNAS2012 import STATES, state

# True: one sigma for all 15 observables (PEtab's Zheng problem: a single estimated
# parameter `sigma` in every noiseFormula).  False: one sigma per observable.
SHARED_SIGMA = True


def no_optimization(replicate):
    return {}


# PEtab's own estimate of sigma (parameters.tsv nominalValue). Levenberg-Marquardt needs every sigma known: the loss is then
# const + 0.5 r.r.  The minimum over the rate constants is the same as with a profiled sigma (it minimizes the same SSE);
# only the value of the objective and the uncertainty estimates differ.
PETAB_SIGMA = 0.002339466852399


def Zheng_loss_config_known_sigma(replicate):
    """As Zheng_loss_config, but with sigma declared instead of estimated (needed by method "LM")."""
    lc = Zheng_loss_config(replicate)
    for o in lc["observables"]:
        o.pop("sigma_pool", None)
        o.update({"sigma_method": "fixed", "sigma_value": PETAB_SIGMA})
    return lc


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

import os

from antimony_modules.Histone.methylation_Zheng_PNAS2012 import STATES, state

# True: one sigma for all 15 observables (PEtab's Zheng problem: a single estimated
# parameter `sigma` in every noiseFormula).  False: one sigma per observable.
SHARED_SIGMA = True


# "column": each observable is a bare result column (the assignment rule observable_K27me<a>K36me<b>).
# "callable": the same number returned by a function of the table. Statistically identical; to the engine it is NOT a bare
# column, which is what a difference of simulations, a normalised sum or a ratio is, so benchmarks/ uses it to time the
# derivative route such observables take.
OBSERVABLE_FORM = os.environ.get("ZHENG_OBSERVABLES", "column")


def _observable(column):
    if OBSERVABLE_FORM == "column":
        return column

    def fraction(result):
        return result[column]
    return fraction


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
        "observed_variable": _observable(f"observable_{state(a, b)}"),
        "data_column":       f"observable_{state(a, b)}",
        "time_column":       "time",
        "data_dict_key":     label,
        # a callable of the table reads no fitted parameter by name: declare it, or every parameter gets a partial
        **({"reads_parameters": []} if OBSERVABLE_FORM != "column" else {}),
        **({"sigma_pool": "sigma"} if SHARED_SIGMA else {}),
    } for a, b in STATES]}

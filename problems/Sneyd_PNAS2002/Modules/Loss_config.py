# PEtab estimates one noise parameter, sigma, shared by every measurement: the loss elements of all nine conditions carry the same
# "sigma_block" (Optimizer_settings.py), so there is one sigma for the 135 points, profiled out of the likelihood.

# PEtab's own value of sigma (parameters.tsv nominalValue). Levenberg-Marquardt needs the sigma known: the loss is then
# const + 0.5 r.r. Its minimum over the rate constants is that of the profiled problem at a self-consistent sigma.
PETAB_SIGMA = 0.022665679544462


def no_optimization(replicate):
    return {}


def Sneyd_loss_config(replicate):
    """PEtab's one observable: the open probability, a rule of the model, with normal noise."""
    return {"observables": [{
        "observed_variable": "open_probability",
        "data_column":       "open_probability",
        "time_column":       "time",
        "data_dict_key":     replicate["Label"],
    }]}


def Sneyd_loss_config_known_sigma(replicate):
    """As Sneyd_loss_config, with sigma declared instead of estimated (needed by method "LM")."""
    lc = Sneyd_loss_config(replicate)
    for o in lc["observables"]:
        o.update({"sigma_method": "fixed", "sigma_value": PETAB_SIGMA})
    return lc

def no_optimization(replicate):
    return {}


def Zhao_loss_config(replicate):
    """PEtab's observable_confirmed_infected: the confirmed-infected species, normal noise.

    Each replicate is its own likelihood block, so each (region, stage) estimates its own sigma: PEtab's per-dataset
    sd_confirmed_infected_<region>_<stage> parameters.
    """
    return {"observables": [{
        "observed_variable": "Confirmed_Infected_Population",
        "data_column":       "Confirmed_Infected",
        "time_column":       "time",
        "data_dict_key":     replicate["Label"],
    }]}

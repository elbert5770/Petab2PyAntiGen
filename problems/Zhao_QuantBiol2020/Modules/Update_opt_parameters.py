# Update_opt_parameters.py: called at every optimization step, to move parameters that are not optimized themselves but
# depend on one that is.  The fitted parameters are named per condition (R_Stage_I_Wuhan, gamma_1_Stage_I_Wuhan, ...);
# the model has one R, gamma_1 and gamma_2, which each replicate takes from its own region and stage.

FITTED = ("R", "gamma_1", "gamma_2")


def update_opt_Zhao(r, replicate, parameters):
    for name in FITTED:
        fitted = f"{name}_{replicate['stage']}_{replicate['region']}"
        if fitted in parameters:
            r[name] = parameters[fitted]

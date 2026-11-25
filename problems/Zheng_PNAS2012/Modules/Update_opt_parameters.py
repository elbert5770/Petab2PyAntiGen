# Update_opt_parameters.py: called at every optimization step, to move parameters that
# are not optimized but depend on one that is.  None here: every parameter that varies
# is itself in the Optimization (Optimizer_settings.py).

def update_opt_no_parameters(r, experiment, parameters):
    return

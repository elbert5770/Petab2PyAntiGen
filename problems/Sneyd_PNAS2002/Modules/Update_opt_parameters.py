# Update_opt_parameters.py: called at every optimization step, to move parameters that are not optimized but depend on one
# that is.  None here: every fitted rate constant is itself a parameter of the model.

def update_opt_no_parameters(r, experiment, parameters):
    return

# Update_opt_parameters.py: called at every optimization step, to move parameters that are not optimized but depend on one
# that is.  None here: every parameter that varies, the observable scale factors included, is itself in the Optimization.

def update_opt_no_parameters(r, experiment, parameters):
    return

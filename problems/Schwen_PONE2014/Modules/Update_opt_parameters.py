# Update_opt_parameters.py: called at every optimization step, to move what depends on a fitted parameter.
# The receptor numbers are initial amounts that PEtab computes from two fitted parameters (Rec1 = ini_R1, Rec2 = ini_R1 * ini_R2fold);
# a species is not a parameter, so this sets it, as it did for Zhao's regional epidemics.

def update_opt_Schwen(r, replicate, parameters):
    ini_R1 = parameters.get("ini_R1", r["ini_R1"])
    ini_R2fold = parameters.get("ini_R2fold", r["ini_R2fold"])
    r["Rec1_ExtracellularMedium"] = ini_R1
    r["Rec2_ExtracellularMedium"] = ini_R1 * ini_R2fold

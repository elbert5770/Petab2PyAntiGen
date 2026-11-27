# Schwen_PONE2014: the fitting routes compared

Every route starts from the same perturbed point (seed 1: each parameter x 10**U(-0.3, 0.3), inside its bounds; a route listed with another seed starts from that seed's point, the same for every route with that seed) and fits the problem's parameters on a log10 scale. Each fitted parameter vector is then **re-scored with one objective, PEtab's own** (SBML + tables, RoadRunner, each noise parameter at its ML value), so the routes are compared on the same footing.

PEtab's published parameters score NLL = -454.4529 (-log L = 956.5184; 2 noise parameters at ML, median sigma = 0.147754). NLL is in the engine's convention, `(n/2) ln(SSE/n)` for each estimated noise group and `n ln(sigma) + SSE/(2 sigma^2)` for a declared sigma; -log L adds N/2 ln 2 pi, and N_est/2 for the points of the estimated groups.

| route | start | sigma | converged | evals | wall s | s/eval | NLL (PEtab reference) | vs PEtab's | median / max shift in decades | > 0.1 decade |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| PEtab -> AMICI directly / L-BFGS-B | seed 1 | estimated jointly (PEtab parameters) | yes | 303 | 108 | 0.358 | -454.4523 | +0.0006 | 0.001 / 0.13 | 2 of 28 |
| PyAntiGen / amici / L-BFGS-B | seed 1 | profiled | yes | 253 | 130 | 0.513 | -454.4513 | +0.0016 | 0.001 / 0.13 | 2 of 28 |
| PyAntiGen / roadrunner / Nelder-Mead | seed 1 | profiled | no (budget) | 20000 | 2894 | 0.145 | -453.9927 | +0.4603 | 0.034 / 3.46 | 9 of 28 |
| PyAntiGen / amici / L-BFGS-B | seed 2 | profiled | yes | 199 | 102 | 0.513 | -454.4475 | +0.0054 | 0.002 / 0.18 | 3 of 28 |
| PEtab -> AMICI directly / L-BFGS-B | seed 2 | estimated jointly (PEtab parameters) | yes | 300 | 106 | 0.354 | -454.4523 | +0.0006 | 0.000 / 0.14 | 1 of 28 |
| PyAntiGen / roadrunner / Nelder-Mead | seed 2 | profiled | no (budget) | 20000 | 2897 | 0.145 | -454.1882 | +0.2648 | 0.018 / 3.62 | 8 of 28 |
| PyAntiGen / amici / L-BFGS-B | seed 3 | profiled | yes | 96 | 52 | 0.546 | -454.4465 | +0.0064 | 0.002 / 0.30 | 3 of 28 |
| PEtab -> AMICI directly / L-BFGS-B | seed 3 | estimated jointly (PEtab parameters) | yes | 190 | 68 | 0.359 | -454.4515 | +0.0014 | 0.002 / 0.30 | 3 of 28 |

Platform / versions: Linux-5.15.167.4-microsoft-standard-WSL2-x86_64-with-glibc2.39, Python 3.12.3

Notes: **direct_amici_lbfgsb**: objective = -llh + priors - N/2 ln(2 pi) - N_est/2 (N_est: points of estimated noise parameters); PEtab's objective priors added (AMICI's llh has none); solver rtol 1e-08, atol 1e-10, retry [(1e-06, 1e-08), (1e-05, 1e-07), (0.0001, 1e-06)] **direct_amici_lbfgsb_seed2**: objective = -llh + priors - N/2 ln(2 pi) - N_est/2 (N_est: points of estimated noise parameters); PEtab's objective priors added (AMICI's llh has none); solver rtol 1e-08, atol 1e-10, retry [(1e-06, 1e-08), (1e-05, 1e-07), (0.0001, 1e-06)] **direct_amici_lbfgsb_seed3**: objective = -llh + priors - N/2 ln(2 pi) - N_est/2 (N_est: points of estimated noise parameters); PEtab's objective priors added (AMICI's llh has none); solver rtol 1e-08, atol 1e-10, retry [(1e-06, 1e-08), (1e-05, 1e-07), (0.0001, 1e-06)]

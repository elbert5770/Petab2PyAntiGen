# Sneyd_PNAS2002: the fitting routes compared

Every route starts from the same perturbed point (seed 1: each parameter x 10**U(-0.3, 0.3), inside its bounds; a route listed with another seed starts from that seed's point, the same for every route with that seed) and fits the problem's parameters on a log10 scale. Each fitted parameter vector is then **re-scored with one objective, PEtab's own** (SBML + tables, RoadRunner, each noise parameter at its ML value), so the routes are compared on the same footing.

PEtab's published parameters score NLL = -511.3486 (-log L = -319.7919; 1 noise parameter at ML, sigma = 0.0226461). NLL is in the engine's convention, the sum over noise groups of `(n/2) ln(SSE/n)`; -log L adds N/2 (1 + ln 2 pi).

| route | start | sigma | converged | evals | wall s | s/eval | NLL (PEtab reference) | vs PEtab's | median / max shift in decades | > 0.1 decade |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| PyAntiGen / amici / LM (chain rule, exact) | seed 1 | declared | yes | 19 | 9 | 0.489 | -511.3490 | -0.0004 | 0.003 / 0.04 | 0 of 14 |
| PEtab -> AMICI directly / L-BFGS-B | seed 1 | estimated jointly (PEtab parameters) | yes | 66 | 19 | 0.289 | -511.3449 | +0.0037 | 0.009 / 0.17 | 1 of 14 |
| PyAntiGen / amici / L-BFGS-B | seed 1 | profiled | yes | 55 | 19 | 0.349 | -511.3449 | +0.0037 | 0.009 / 0.17 | 1 of 14 |
| PyAntiGen / roadrunner / Nelder-Mead | seed 1 | profiled | yes | 4455 | 41 | 0.009 | -511.3236 | +0.0250 | 0.021 / 2.21 | 5 of 14 |
| PyAntiGen / amici / LM (chain rule, exact) | seed 2 | declared | yes | 38 | 14 | 0.366 | -511.3489 | -0.0003 | 0.004 / 0.16 | 2 of 14 |
| PyAntiGen / amici / L-BFGS-B | seed 2 | profiled | yes | 79 | 24 | 0.306 | -511.3488 | -0.0002 | 0.038 / 0.13 | 4 of 14 |
| PEtab -> AMICI directly / L-BFGS-B | seed 2 | estimated jointly (PEtab parameters) | yes | 91 | 26 | 0.286 | -511.3487 | -0.0001 | 0.036 / 0.12 | 4 of 14 |

Platform / versions: Linux-5.15.167.4-microsoft-standard-WSL2-x86_64-with-glibc2.39, Python 3.12.3

Notes: **direct_amici_lbfgsb**: objective = -llh - N/2 ln(2 pi) - N_est/2 (N_est: points of estimated noise parameters) **direct_amici_lbfgsb_seed2**: objective = -llh - N/2 ln(2 pi) - N_est/2 (N_est: points of estimated noise parameters)

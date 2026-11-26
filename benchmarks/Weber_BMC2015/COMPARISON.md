# Weber_BMC2015: the fitting routes compared

Every route starts from the same perturbed point (seed 1: each parameter x 10**U(-0.3, 0.3), inside its bounds; a route listed with another seed starts from that seed's point, the same for every route with that seed) and fits the problem's parameters on a log10 scale. Each fitted parameter vector is then **re-scored with one objective, PEtab's own** (SBML + tables, RoadRunner, each noise parameter at its ML value), so the routes are compared on the same footing.

PEtab's published parameters score NLL = 112.1451 (-log L = 296.2018; 8 noise parameters at ML, median sigma = 1.0272). NLL is in the engine's convention, the sum over noise groups of `(n/2) ln(SSE/n)`; -log L adds N/2 (1 + ln 2 pi).

| route | start | sigma | converged | evals | wall s | s/eval | NLL (PEtab reference) | vs PEtab's | median / max shift in decades | > 0.1 decade |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| PyAntiGen / amici / L-BFGS-B | seed 1 | profiled, some declared | yes | 19 | 8 | 0.437 | reference solver fails (engine: 12485.1) | | 5.628 / 9.34 | 29 of 31 |
| PEtab -> AMICI directly / L-BFGS-B | seed 1 | estimated jointly (PEtab parameters) | yes | 1709 | 338 | 0.198 | 112.1580 | +0.0129 | 0.025 / 0.63 | 10 of 31 |
| PyAntiGen / roadrunner / Nelder-Mead | seed 1 | profiled, some declared | no (budget) | 20000 | 546 | 0.027 | 114.6700 | +2.5249 | 0.118 / 2.77 | 21 of 31 |
| PyAntiGen / amici / L-BFGS-B | seed 2 | profiled, some declared | yes | 312 | 72 | 0.230 | 112.2457 | +0.1006 | 0.045 / 0.24 | 8 of 31 |
| PEtab -> AMICI directly / L-BFGS-B | seed 2 | estimated jointly (PEtab parameters) | yes | 1148 | 235 | 0.204 | 112.1513 | +0.0062 | 0.019 / 0.14 | 1 of 31 |

Platform / versions: Linux-5.15.167.4-microsoft-standard-WSL2-x86_64-with-glibc2.39, Python 3.12.3

Notes: **pyantigen_amici_lbfgsb**: L-BFGS-B reports convergence, but at a corner of the box: all 31 parameters end within 0.5 decade of a bound (20 within 0.05), where the objective is flat. Its first line search stepped 30 in log10 space (the start's gradient norm is 4e5) and dropped the objective from 21,646 to 12,573. The gradient is exact at the start (31 of 31 agree with finite differences of PEtab's reference, max relative error 7e-8), so the optimizer's path is at fault, not the gradient. The reference solver cannot integrate this vector (stiff), so it is not rescored. **direct_amici_lbfgsb**: objective = -llh - N/2 ln(2 pi) - N_est/2 (N_est: points of estimated noise parameters) **direct_amici_lbfgsb_seed2**: objective = -llh - N/2 ln(2 pi) - N_est/2 (N_est: points of estimated noise parameters)

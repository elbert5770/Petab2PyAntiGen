# Zheng_PNAS2012: the fitting routes compared

Every route starts from the same perturbed point (seed 1: each parameter x 10**U(-0.3, 0.3)) and fits all 45 rate constants and `inflowp` on a log10 scale in [1e-5, 1e3]. Each fitted parameter vector is then **re-scored with one objective, PEtab's own** (SBML + tables, RoadRunner, steady-state pre-equilibration, sigma at its ML value), so the routes are compared on the same footing.

PEtab's published parameters score NLL = -363.4698 (-log L = -278.3335, sigma = 0.002339). NLL is in the engine's convention, `(n/2) ln(SSE/n)`; -log L adds n/2 (1 + ln 2 pi).

| route | sigma | converged | evals | wall s | s/eval | NLL (PEtab reference) | vs PEtab's | median / max shift in decades | > 0.1 decade |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| PyAntiGen / amici / LM | declared | yes | 368 | 48 | 0.130 | -363.4692 | +0.0006 | 0.001 / 0.19 | 1 of 45 |
| PEtab -> AMICI directly / L-BFGS-B | estimated jointly (PEtab parameter) | yes | 474 | 60 | 0.126 | -363.4649 | +0.0049 | 0.004 / 0.31 | 11 of 45 |
| PyAntiGen / amici / L-BFGS-B | profiled | yes | 465 | 68 | 0.147 | -363.4650 | +0.0049 | 0.004 / 0.31 | 11 of 45 |
| PyAntiGen / roadrunner / Nelder-Mead | profiled | no (budget) | 20000 | 119 | 0.006 | -363.2856 | +0.1842 | 0.053 / 3.71 | 14 of 45 |

Platform / versions: Linux-5.15.167.4-microsoft-standard-WSL2-x86_64-with-glibc2.39, Python 3.12.3

Notes: **direct_amici_lbfgsb**: objective = -llh - n/2 (1 + ln 2 pi)

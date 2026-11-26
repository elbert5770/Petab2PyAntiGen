# Zhao_QuantBiol2020: the fitting routes compared

Every route starts from the same perturbed point (seed 1: each parameter x 10**U(-0.3, 0.3), inside its bounds) and fits the problem's parameters on a log10 scale. Each fitted parameter vector is then **re-scored with one objective, PEtab's own** (SBML + tables, RoadRunner, each noise parameter at its ML value), so the routes are compared on the same footing.

PEtab's published parameters score NLL = 520.6265 (-log L = 636.9795; 7 noise parameters at ML, median sigma = 481.907). NLL is in the engine's convention, the sum over noise groups of `(n/2) ln(SSE/n)`; -log L adds N/2 (1 + ln 2 pi).

| route | sigma | converged | evals | wall s | s/eval | NLL (PEtab reference) | vs PEtab's | median / max shift in decades | > 0.1 decade |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| PyAntiGen / roadrunner / Nelder-Mead | profiled | yes | 6455 | 63 | 0.010 | 426.4508 | -94.1757 | 0.050 / 1.01 | 7 of 21 |
| PEtab -> AMICI directly / L-BFGS-B | estimated jointly (PEtab parameters) | yes | 825 | 135 | 0.164 | 390.7382 | -129.8883 | 0.088 / 2.29 | 10 of 21 |
| PyAntiGen / amici / L-BFGS-B | profiled | yes | 1346 | 197 | 0.146 | 387.5194 | -133.1072 | 0.088 / 3.20 | 10 of 21 |

Platform / versions: Linux-5.15.167.4-microsoft-standard-WSL2-x86_64-with-glibc2.39, Python 3.12.3

Notes: **direct_amici_lbfgsb**: objective = -llh - n/2 (1 + ln 2 pi)

# Smith_BMCSystBiol2013: the fitting routes compared

Every route starts from the same perturbed point (seed 1: each parameter x 10**U(-0.3, 0.3), inside its bounds; a route listed with another seed starts from that seed's point, the same for every route with that seed) and fits the problem's parameters on a log10 scale. Each fitted parameter vector is then **re-scored with one objective, PEtab's own** (SBML + tables, RoadRunner, each noise parameter at its ML value), so the routes are compared on the same footing.

PEtab's published parameters score NLL = 888141.4648 (-log L = 888198.4390; 1 noise parameter at ML, sigma = 1). NLL is in the engine's convention, `(n/2) ln(SSE/n)` for each estimated noise group and `n ln(sigma) + SSE/(2 sigma^2)` for a declared sigma; -log L adds N/2 ln 2 pi, and N_est/2 for the points of the estimated groups.

| route | start | sigma | converged | evals | wall s | s/eval | NLL (PEtab reference) | vs PEtab's | median / max shift in decades | > 0.1 decade |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| PyAntiGen / amici / L-BFGS-B | seed 1 | declared | no (budget) | 13 | 117 | 9.000 | 639477.5421 | -248663.9227 | 0.143 / 0.29 | 17 of 25 |
| PyAntiGen / amici / LM (chain rule, exact) | seed 1 | declared | yes | 107 | 668 | 6.241 | 24535.3698 | -863606.0950 | 1.128 / 2.00 | 20 of 25 |
| PyAntiGen / amici / LM / smoothed switch (chain rule, exact) | seed 1 | declared | yes | 151 | 1564 | 10.355 | 24533.1158 | -863608.3490 | 1.122 / 2.00 | 20 of 25 |
| PyAntiGen / roadrunner / Nelder-Mead | seed 1 | declared | no (budget) | 3000 | 3468 | 1.156 | 27782.1102 | -860359.3545 | 0.352 / 1.99 | 20 of 25 |
| PEtab -> AMICI directly / L-BFGS-B | seed 1 | estimated jointly (PEtab parameters) | yes | 1002 | 5706 | 5.695 | 24535.5561 | -863605.9087 | 0.493 / 2.00 | 19 of 25 |
| PyAntiGen / amici / L-BFGS-B / smoothed switch | seed 1 | declared | yes | 817 | 6947 | 8.503 | 24536.8657 | -863604.5991 | 0.543 / 2.00 | 19 of 25 |
| PyAntiGen / amici / L-BFGS-B | seed 2 | declared | yes | 478 | 1034 | 2.162 | 28867.9713 | -859273.4934 | 0.117 / 1.26 | 13 of 25 |
| PyAntiGen / amici / LM (chain rule, exact) | seed 2 | declared | no (budget) | 182 | 1151 | 6.324 | reference solver fails (engine: 23542.0) | | 0.455 / 2.00 | 21 of 25 |
| PEtab -> AMICI directly / L-BFGS-B | seed 2 | estimated jointly (PEtab parameters) | yes | 279 | 1675 | 6.003 | 24985.0248 | -863156.4399 | 0.182 / 2.00 | 17 of 25 |
| PyAntiGen / amici / LM / smoothed switch (chain rule, exact) | seed 2 | declared | yes | 183 | 1919 | 10.487 | 20865.1894 | -867276.2753 | 0.674 / 2.00 | 21 of 25 |
| PyAntiGen / amici / L-BFGS-B / smoothed switch | seed 2 | declared | yes | 335 | 3353 | 10.010 | 24954.2085 | -863187.2563 | 0.212 / 2.00 | 17 of 25 |

Platform / versions: Linux-5.15.167.4-microsoft-standard-WSL2-x86_64-with-glibc2.39, Python 3.12.3; Windows-11-10.0.26200-SP0, Python 3.13.3

Notes: **pyantigen_amici_lm_smooth**: smoothed Akt/PKC switch, eps 0.1 molecules (not PEtab's model); scored on the exact PEtab model **direct_amici_lbfgsb**: objective = -llh - N/2 ln(2 pi) - N_est/2 (N_est: points of estimated noise parameters); solver rtol 1e-08, atol 1e-10, retry [(1e-06, 1e-08), (1e-05, 1e-07), (0.0001, 1e-06)] **pyantigen_amici_lbfgsb_smooth**: smoothed Akt/PKC switch, eps 0.1 molecules (not PEtab's model); scored on the exact PEtab model **pyantigen_amici_lm_seed2**: the engine discards a fit it judges not converged; these are its best parameters so far, from its progress snapshot (results/*optimization_progress_20261126_052321.json); wall time estimated from the log's timestamps (the run started together with the L-BFGS-B start 2 fit); AMICI's own likelihood (PEtab importer, rtol 1e-10) at these parameters gives NLL 23542.04, the engine's value, but the RoadRunner reference cannot integrate them **direct_amici_lbfgsb_seed2**: objective = -llh - N/2 ln(2 pi) - N_est/2 (N_est: points of estimated noise parameters); solver rtol 1e-08, atol 1e-10, retry [(1e-06, 1e-08), (1e-05, 1e-07), (0.0001, 1e-06)] **pyantigen_amici_lm_smooth_seed2**: smoothed Akt/PKC switch, eps 0.1 molecules (not PEtab's model); scored on the exact PEtab model **pyantigen_amici_lbfgsb_smooth_seed2**: smoothed Akt/PKC switch, eps 0.1 molecules (not PEtab's model); scored on the exact PEtab model

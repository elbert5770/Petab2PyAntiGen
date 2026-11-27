#!/bin/bash
# Fit Schwen_PONE2014 by every route, then compare them (WSL / Linux with AMICI).
#   PYANTIGEN_DIR  the PyAntiGen checkout (default: a sibling of this repository)
#   PY             Python with AMICI and PyAntiGen's requirements   (default /opt/venv/bin/python)
#   PY_PETAB       the same plus the `petab` package, for the direct route (default /root/venv_petab/bin/python)
# Three perturbed starts (seeds 1-3) for the AMICI routes, two (seeds 1 and 2) for Nelder-Mead (20,000 evaluations take ~55 min each).
# No Levenberg-Marquardt route: PEtab estimates two noise parameters, and
# least squares needs every sigma declared. Nelder-Mead on RoadRunner has no priors (they are part of the AMICI route's objective);
# compare_fits.py scores every fit with PEtab's objective, priors included.
set -e
HERE="$(cd "$(dirname "$0")/../.." && pwd)"
export PYANTIGEN_DIR="${PYANTIGEN_DIR:-$(dirname "$HERE")/PyAntiGen}"
export PYTHONPATH="$PYANTIGEN_DIR"
export PATH="$(dirname "${PY:-/opt/venv/bin/python}"):$PATH"       # AMICI's own SWIG must come first
PY="${PY:-/opt/venv/bin/python}"
PY_PETAB="${PY_PETAB:-/root/venv_petab/bin/python}"
cd "$HERE"
P=Schwen_PONE2014

for SEED in 1 2 3; do
  SUF=""; [ "$SEED" != "1" ] && SUF="_seed$SEED"
  "$PY" -u benchmarks/fit_problem.py --problem $P --label pyantigen_amici_lbfgsb$SUF --spec OPTIMIZATION_Schwen_amici --seed $SEED
  "$PY_PETAB" -u benchmarks/amici_direct.py --problem $P --spec OPTIMIZATION_Schwen --seed $SEED --rtol 1e-8 --atol 1e-10 --retry 1e-6:1e-8,1e-5:1e-7,1e-4:1e-6
done
for SEED in 1 2; do
  SUF=""; [ "$SEED" != "1" ] && SUF="_seed$SEED"
  "$PY" -u benchmarks/fit_problem.py --problem $P --label pyantigen_roadrunner_neldermead$SUF --spec OPTIMIZATION_Schwen --seed $SEED --max-evals 20000
done
echo "all fits done; compare them (on any platform with the tool's venv):  python benchmarks/compare_fits.py --problem $P"

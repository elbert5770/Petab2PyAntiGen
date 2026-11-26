#!/bin/bash
# Fit Sneyd_PNAS2002 by every route, then compare them (WSL / Linux with AMICI).
#   PYANTIGEN_DIR  the PyAntiGen checkout (default: a sibling of this repository)
#   PY             Python with AMICI and PyAntiGen's requirements   (default /opt/venv/bin/python)
#   PY_PETAB       the same plus the `petab` package, for the direct route (default /root/venv_petab/bin/python)
# The AMICI routes run from two perturbed starts (seeds 1 and 2): on Weber one start said little about a route.
set -e
HERE="$(cd "$(dirname "$0")/../.." && pwd)"
export PYANTIGEN_DIR="${PYANTIGEN_DIR:-$(dirname "$HERE")/PyAntiGen}"
export PYTHONPATH="$PYANTIGEN_DIR"
export PATH="$(dirname "${PY:-/opt/venv/bin/python}"):$PATH"       # AMICI's own SWIG must come first
PY="${PY:-/opt/venv/bin/python}"
PY_PETAB="${PY_PETAB:-/root/venv_petab/bin/python}"
cd "$HERE"
P=Sneyd_PNAS2002

"$PY" -u benchmarks/fit_problem.py --problem $P --label pyantigen_roadrunner_neldermead --spec OPTIMIZATION_Sneyd
for SEED in 1 2; do
  SUF=""; [ "$SEED" != "1" ] && SUF="_seed$SEED"
  "$PY" -u benchmarks/fit_problem.py --problem $P --label pyantigen_amici_lbfgsb$SUF --spec OPTIMIZATION_Sneyd_amici --seed $SEED
  "$PY" -u benchmarks/fit_problem.py --problem $P --label pyantigen_amici_lm$SUF     --spec OPTIMIZATION_Sneyd_lm    --seed $SEED
  "$PY_PETAB" -u benchmarks/amici_direct.py --problem $P --spec OPTIMIZATION_Sneyd --seed $SEED
done
echo "all fits done; compare them (on any platform with the tool's venv):  python benchmarks/compare_fits.py --problem $P"

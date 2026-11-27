#!/bin/bash
# Fit Smith_BMCSystBiol2013 by every route, then compare them (WSL / Linux with AMICI).
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
P=Smith_BMCSystBiol2013

"$PY" -u benchmarks/fit_problem.py --problem $P --label pyantigen_roadrunner_neldermead --spec OPTIMIZATION_Smith --max-evals 3000   # ~1.4 s per evaluation: 20,000 would take hours
for SEED in 1 2; do
  SUF=""; [ "$SEED" != "1" ] && SUF="_seed$SEED"
  "$PY" -u benchmarks/fit_problem.py --problem $P --label pyantigen_amici_lbfgsb$SUF --spec OPTIMIZATION_Smith_amici --seed $SEED
  "$PY" -u benchmarks/fit_problem.py --problem $P --label pyantigen_amici_lm$SUF     --spec OPTIMIZATION_Smith_lm    --seed $SEED
  "$PY_PETAB" -u benchmarks/amici_direct.py --problem $P --spec OPTIMIZATION_Smith --seed $SEED --rtol 1e-8 --atol 1e-10 --retry 1e-6:1e-8,1e-5:1e-7,1e-4:1e-6 --max-evals 1500
done
# the smoothed Akt/PKC switch (an AMICI-friendly variant of the model; the fits are scored on PEtab's exact model by compare_fits.py).
# Run it after the exact fits above: it regenerates the model files, so leave nothing else running, and regenerate the exact files after.
export SMITH_PIP3_SWITCH_EPS=0.1
for SEED in 1 2; do
  SUF=""; [ "$SEED" != "1" ] && SUF="_seed$SEED"
  "$PY" -u benchmarks/fit_problem.py --problem $P --label pyantigen_amici_lbfgsb_smooth$SUF --spec OPTIMIZATION_Smith_amici --seed $SEED --note "smoothed Akt/PKC switch, eps 0.1 molecules"
  "$PY" -u benchmarks/fit_problem.py --problem $P --label pyantigen_amici_lm_smooth$SUF     --spec OPTIMIZATION_Smith_lm    --seed $SEED --note "smoothed Akt/PKC switch, eps 0.1 molecules"
done
unset SMITH_PIP3_SWITCH_EPS
(cd problems/$P && "$PY" Model_generate.py > /dev/null)         # the exact model again
echo "all fits done; compare them (on any platform with the tool's venv):  python benchmarks/compare_fits.py --problem $P"

#!/bin/bash
# Fit Zheng_PNAS2012 by every route, from the same perturbed start, then compare them (WSL / Linux with AMICI).
#   PYANTIGEN_DIR  the PyAntiGen checkout (default: a sibling of this repository)
#   PY             Python with AMICI and PyAntiGen's requirements   (default /opt/venv/bin/python)
#   PY_PETAB       the same plus the `petab` package, for the direct route (default /root/venv_petab/bin/python)
set -e
HERE="$(cd "$(dirname "$0")/.." && pwd)"
export PYANTIGEN_DIR="${PYANTIGEN_DIR:-$(dirname "$HERE")/PyAntiGen}"
export PYTHONPATH="$PYANTIGEN_DIR"
export PATH="$(dirname "${PY:-/opt/venv/bin/python}"):$PATH"       # AMICI's own SWIG must come first
PY="${PY:-/opt/venv/bin/python}"
PY_PETAB="${PY_PETAB:-/root/venv_petab/bin/python}"
cd "$HERE"

"$PY" -u benchmarks/fit_zheng.py --label pyantigen_roadrunner_neldermead --spec OPTIMIZATION_Zheng
"$PY" -u benchmarks/fit_zheng.py --label pyantigen_amici_lbfgsb          --spec OPTIMIZATION_Zheng_amici
"$PY" -u benchmarks/fit_zheng.py --label pyantigen_amici_lm              --spec OPTIMIZATION_Zheng_lm
# the same LM fit with finite-difference Jacobians: what an observable that is not a bare column gets (about 5 minutes)
"$PY" -u benchmarks/fit_zheng.py --label pyantigen_amici_lm_finite_difference --spec OPTIMIZATION_Zheng_lm --no-exact
# the same LM fit with observables that are callables of the table, not bare columns (as a difference of simulations is):
# the chain rule differences only the observable; the old route differenced the whole loss
"$PY" -u benchmarks/fit_zheng.py --label pyantigen_amici_lm_callable_observable --spec OPTIMIZATION_Zheng_lm --observables callable
"$PY" -u benchmarks/fit_zheng.py --label pyantigen_amici_lm_callable_finite_difference --spec OPTIMIZATION_Zheng_lm --observables callable --no-exact
"$PY_PETAB" -u benchmarks/amici_direct_zheng.py
echo "all fits done; compare them (on any platform with the tool's venv):  python benchmarks/compare_fits.py"

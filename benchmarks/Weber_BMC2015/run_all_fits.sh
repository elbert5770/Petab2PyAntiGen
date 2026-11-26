#!/bin/bash
# Fit Weber_BMC2015 by every route, from the same perturbed start, then compare them (WSL / Linux with AMICI).
#   PYANTIGEN_DIR  the PyAntiGen checkout (default: a sibling of this repository)
#   PY             Python with AMICI and PyAntiGen's requirements   (default /opt/venv/bin/python)
#   PY_PETAB       the same plus the `petab` package, for the direct route (default /root/venv_petab/bin/python)
set -e
HERE="$(cd "$(dirname "$0")/../.." && pwd)"
export PYANTIGEN_DIR="${PYANTIGEN_DIR:-$(dirname "$HERE")/PyAntiGen}"
export PYTHONPATH="$PYANTIGEN_DIR"
export PATH="$(dirname "${PY:-/opt/venv/bin/python}"):$PATH"       # AMICI's own SWIG must come first
PY="${PY:-/opt/venv/bin/python}"
PY_PETAB="${PY_PETAB:-/root/venv_petab/bin/python}"
cd "$HERE"

"$PY" -u benchmarks/fit_problem.py --problem Weber_BMC2015 --label pyantigen_roadrunner_neldermead --spec OPTIMIZATION_Weber
"$PY" -u benchmarks/fit_problem.py --problem Weber_BMC2015 --label pyantigen_amici_lbfgsb          --spec OPTIMIZATION_Weber_amici
"$PY_PETAB" -u benchmarks/amici_direct.py --problem Weber_BMC2015 --spec OPTIMIZATION_Weber
# A second perturbed start for the two AMICI routes: from seed 1 (the start every route shares) L-BFGS-B's first line search steps
# to a corner of the box on PyAntiGen's objective and stops there, so one start says little about either route.
"$PY" -u benchmarks/fit_problem.py --problem Weber_BMC2015 --label pyantigen_amici_lbfgsb_seed2 --spec OPTIMIZATION_Weber_amici --seed 2
"$PY_PETAB" -u benchmarks/amici_direct.py --problem Weber_BMC2015 --spec OPTIMIZATION_Weber --seed 2
echo "all fits done; compare them (on any platform with the tool's venv):  python benchmarks/compare_fits.py --problem Weber_BMC2015"

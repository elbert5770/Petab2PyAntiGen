#!/bin/bash
# Fit Zhao_QuantBiol2020 by every route, from the same perturbed start, then compare them (WSL / Linux with AMICI).
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

"$PY" -u benchmarks/fit_problem.py --problem Zhao_QuantBiol2020 --label pyantigen_roadrunner_neldermead --spec OPTIMIZATION_Zhao
"$PY" -u benchmarks/fit_problem.py --problem Zhao_QuantBiol2020 --label pyantigen_amici_lbfgsb          --spec OPTIMIZATION_Zhao_amici
"$PY_PETAB" -u benchmarks/amici_direct.py --problem Zhao_QuantBiol2020 --spec OPTIMIZATION_Zhao
echo "all fits done; compare them (on any platform with the tool's venv):  python benchmarks/compare_fits.py --problem Zhao_QuantBiol2020"

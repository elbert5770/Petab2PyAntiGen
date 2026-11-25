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
"$PY_PETAB" -u benchmarks/amici_direct_zheng.py
echo "all fits done; compare them (on any platform with the tool's venv):  python benchmarks/compare_fits.py"

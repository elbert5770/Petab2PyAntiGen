# Petab2PyAntiGen

[PEtab](https://petab.readthedocs.io) states a parameter-estimation problem as an SBML file and
five or six TSV tables. For regular models that is huge: the SBML of the benchmark problem
Zheng_PNAS2012 is 95 kB for what is, in the paper, a few rules.
[PyAntiGen](https://github.com/elbert5770/PyAntiGen) writes the same problem as Python
f-strings that loop over species, isotopes and compartments, keeps compartment volumes explicit,
names species `<species>_<Compartment>`, and puts the estimation problem in
`Optimizer_settings.py` / `Loss_config.py` / `Update_opt_parameters.py` / `Update_parameters.py`.

This repository is a tool for doing that conversion, plus the converted problems.

```
petab2pyantigen/    the tool: PEtab reader, reference simulator, SBML->Antimony, emit, verify, stats
problems/<Name>/    one PyAntiGen project per converted problem (Model_run.py, Modules/, ...)
antimony_modules/   the reaction modules (f-string loops)
antimony_models/    parameters, initial conditions, manual rules, events  (+ generated reactions)
data/               measurements, wide format
tests/              each converted problem is tested against PEtab
.agents/skills/     how to convert a problem (worked example: Zheng_PNAS2012)
```

## Setup

Python >= 3.11. Everything is installed into `.venv/` (git-ignored).

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python -m pip install --no-deps -e .
```

## The tool

```powershell
# SBML -> Antimony through Tellurium, minus `X is "..."` lines, annotations and notes
python -m petab2pyantigen antimony <PEtab dir>\model.xml

# parameters.csv, InitialConditions.csv and the wide data csv from a PEtab problem
python -m petab2pyantigen emit <PEtab dir> --model Zheng_PNAS2012

# does the PyAntiGen project reproduce PEtab? (nominal + random parameters)
python -m petab2pyantigen verify --petab <PEtab dir> --project problems/Zheng_PNAS2012 --experiment Zheng

# size of the PEtab problem vs the PyAntiGen project
python -m petab2pyantigen stats <PEtab dir> --model Zheng_PNAS2012 --modules Histone/methylation_Zheng_PNAS2012.py
```

The PEtab problems themselves are not copied here; point at a checkout of the PEtab benchmark
collection (tests look for `../petab-benchmark-problems`, or `PETAB_BENCHMARK_DIR`).

## Zheng_PNAS2012

Combinatorial methylation of histone H3 K27 and K36 (Zheng et al., PNAS 2012): 15 states on a
4 x 4 lattice, 44 rate constants, inflow/dilution, a pre-equilibration then a 25 h time course.

| | PEtab | PyAntiGen |
| --- | --- | --- |
| reactions, rules | 2,832-line SBML (95 kB), 60 reaction function definitions | `antimony_modules/Histone/methylation_Zheng_PNAS2012.py`, 3 loops |
| observables, noise | `observables.tsv`, 15 long formulas (5 kB) | one loop in `Loss_config.py` + assignment rules |
| conditions | `conditions.tsv` + pre-equilibration | one event, one untracked block |
| data | 61-row long table | 4 rows x 15 columns |

```powershell
cd problems\Zheng_PNAS2012
python Model_run.py --simulate Zheng              # figure in results/Zheng_PNAS2012/
python Model_run.py --optimize Zheng --no-fit     # objective at PEtab's nominal parameters
python Model_run.py --optimize Zheng              # fit all 45 parameters (log10, [1e-5, 1e3])
```

Setup note: the single shared sigma uses PyAntiGen's `"sigma_pool"` option, which is newer than
the 2.1.6 release; install the PyAntiGen checkout with `pip install --no-deps -e ../PyAntiGen`.

Checked against PEtab (`tests/test_zheng.py`): observables agree to ~1e-11 at the nominal and at
random parameters, and the objective at the nominal parameters is PEtab's: one pooled sigma =
0.0023395, -log L = -278.3335 (the engine prints -363.47, the same number without the
constant n/2 (1 + ln 2 pi)).

### AMICI (WSL)

AMICI does not run on native Windows. In WSL (the default `ubuntu-amici` distro has `/opt/venv` with AMICI 1.1.0
and the toolchain) run the same project with the venv first on `PATH`, so that AMICI's own SWIG (4.5.1) is used,
not the system one (4.2.0, which compiles a model that then fails to load), and the PyAntiGen checkout on `PYTHONPATH`:

```bash
export PATH=/opt/venv/bin:$PATH PYTHONPATH=/mnt/c/Users/elber/Documents/git/PyAntiGen
cd /mnt/c/Users/elber/Documents/git/Petab2PyAntiGen/problems/Zheng_PNAS2012
python Model_run.py --simulate Zheng --backend amici --compare   # AMICI vs RoadRunner: 1e-7 relative
python Model_run.py --optimize ZhengAMICI                        # L-BFGS-B with the analytic gradient
python Model_run.py --optimize ZhengLM                           # Levenberg-Marquardt (sigma declared: LM needs known sigmas)
```

The model compiles in about a minute (cached in `~/.cache/pyantigen/amici_models`). The AMICI trajectories
agree with RoadRunner to 1.1e-7 relative, including the 1e5 h pre-equilibration and the event.

### Comparing the routes

`benchmarks/` fits Zheng by every route from the same perturbed start (seed 1) and scores each fitted parameter vector
with one objective, PEtab's own (SBML + tables, RoadRunner, steady-state pre-equilibration), so the routes are compared on the
same footing whatever engine produced them:

```bash
benchmarks/run_all_fits.sh                 # WSL/Linux: the four fits (about 5 minutes in total), writes benchmarks/fits/*.json
python benchmarks/compare_fits.py          # any platform: writes benchmarks/COMPARISON.md and comparison.png
```

* **`benchmarks/COMPARISON.md`** is the table (evaluations, wall-clock, NLL under the PEtab reference, how far the
  parameters moved from PEtab's published ones); **`benchmarks/comparison.png`** plots time against fit quality.
* **`benchmarks/fits/<route>.json`** holds each route's fitted parameters, so a route can be re-scored or compared pairwise.
* The routes: PyAntiGen on RoadRunner (Nelder-Mead), PyAntiGen on AMICI (L-BFGS-B with the gradient; Levenberg-Marquardt
  with sigma declared), and PEtab straight through AMICI with no PyAntiGen (`amici_direct_zheng.py`).
* Single-model checks: `Model_run.py --simulate Zheng --backend amici --compare` prints AMICI against RoadRunner column by
  column (1.1e-7 relative); `python -m petab2pyantigen verify ...` compares PyAntiGen with the PEtab reference.

PyAntiGen's AMICI gradient used to be assembled from ~580 finite-difference re-evaluations of its Python loss (1.9 s) and a
15-panel progress figure was redrawn at every improvement; AMICI itself was under 2% of the time. PyAntiGen now takes the
exact chain-rule gradient (and Levenberg-Marquardt Jacobian) for observables that are bare result columns, and draws the
progress figure once, when the fit ends (`PYANTIGEN_PROGRESS_FIGURE=live` brings back the refreshing one).

Rate constants are named as in PEtab (`k00_01` = K27me0K36me0 -> K27me0K36me1); the 15 observables
are `observable_K27me<a>K36me<b>`.

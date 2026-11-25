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

Rate constants are named as in PEtab (`k00_01` = K27me0K36me0 -> K27me0K36me1); the 15 observables
are `observable_K27me<a>K36me<b>`.

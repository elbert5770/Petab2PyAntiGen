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

### How much does it compress?

`python benchmarks/Zheng_PNAS2012/compression.py` writes `benchmarks/Zheng_PNAS2012/COMPRESSION.md` and `compression.png`: the model (SBML, Tellurium's
Antimony, the Antimony PyAntiGen generates, the module that generates it) and, separately, the estimation problem (each PEtab
table against the PyAntiGen files that replace it), by bytes, code lines and gzip, with the module's growth from a 2 x 2 to a
12 x 12 lattice. In short: the module is 31x smaller than the SBML (8x after gzip) and does not grow with the system; the
generated Antimony is no shorter than Tellurium's; the estimation problem is about the same size in bytes as PEtab's tables
(48 code lines state it, against 156 table lines) and larger for conditions.

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

`benchmarks/Zheng_PNAS2012/` fits Zheng by every route from the same perturbed start (seed 1) and scores each fitted parameter vector
with one objective, PEtab's own (SBML + tables, RoadRunner, steady-state pre-equilibration), so the routes are compared on the
same footing whatever engine produced them:

```bash
benchmarks/Zheng_PNAS2012/run_all_fits.sh   # WSL/Linux: the fits (about 15 minutes in total), writes benchmarks/Zheng_PNAS2012/fits/*.json
python benchmarks/compare_fits.py --problem Zheng_PNAS2012   # any platform: writes COMPARISON.md and comparison.png there
python benchmarks/plot_fits.py --problem Zheng_PNAS2012      # writes fits.png there: every route's final fit over the data
```

* **`benchmarks/Zheng_PNAS2012/COMPARISON.md`** is the table (evaluations, wall-clock, NLL under the PEtab reference, how far the
  parameters moved from PEtab's published ones); **`comparison.png`** plots time to fit by route.
* **`fits/<route>.json`** (next to it) holds each route's fitted parameters, so a route can be re-scored or compared pairwise.
* The routes: PyAntiGen on RoadRunner (Nelder-Mead), PyAntiGen on AMICI (L-BFGS-B with the gradient; Levenberg-Marquardt with
  sigma declared), and PEtab straight through AMICI with no PyAntiGen (`benchmarks/amici_direct.py`). The LM fit is run three ways,
  because the derivative route depends on the observables: **bare result columns** (Zheng's own: exact chain rule, 48 s),
  **callables of the table** (not bare, as a difference of simulations is: only the observable is differenced, 67 s) and
  **finite differences of the whole loss** (what non-bare observables took before, 294-303 s). All three converge to the same
  NLL. Zheng has 45 parameters; a model with 9-18 has proportionally fewer loss passes in the old route, so its saving is smaller.
* Single-model checks: `Model_run.py --simulate Zheng --backend amici --compare` prints AMICI against RoadRunner column by
  column (1.1e-7 relative); `python -m petab2pyantigen verify ...` compares PyAntiGen with the PEtab reference.

PyAntiGen's AMICI gradient used to be assembled from ~580 finite-difference re-evaluations of its Python loss (1.9 s) and a
15-panel progress figure was redrawn at every improvement; AMICI itself was under 2% of the time. PyAntiGen now takes the
chain-rule gradient (and Levenberg-Marquardt Jacobian): exact for observables that are bare result columns, and with only the
observable (or, for a composite, its aggregation) differenced otherwise; it draws the progress figure once, when the fit ends (`PYANTIGEN_PROGRESS_FIGURE=live` brings back the refreshing one).

Rate constants are named as in PEtab (`k00_01` = K27me0K36me0 -> K27me0K36me1); the 15 observables
are `observable_K27me<a>K36me<b>`.

## Zhao_QuantBiol2020

COVID-19 in Wuhan, Hubei and the rest of China (Zhao et al., Quantitative Biology 2020): a three-reaction chain
`Susceptible -> Unquarantined -> Quarantined -> Confirmed`, fitted to the cumulative confirmed cases of seven datasets
(three regions x up to three stages of control). It is the opposite of Zheng: the model is tiny and the weight is in the
**conditions**. PEtab selects each dataset's rates, population and initial infections with six 0/1 trigger parameters
multiplied into SBML rules (`R := Trigger_Wuhan * (Trigger_Stage_I * R_Stage_I_Wuhan + ...)`, and the same for gamma_1,
gamma_2, the population and four initial assignments). Here:

* `antimony_modules/Epidemic/siqc_Zhao_QuantBiol2020.py`: the three reactions (one `custom_conc_per_time` law, two `MA`) and
  the rules `alpha := R * gamma_1`, `beta := ...`.
* `data/Zhao_QuantBiol2020_conditions.csv`: one row per (region, stage): population, initial infections and PEtab's nominal
  `R, gamma_1, gamma_2`. **The tool derives it:** `emit.effective_condition_values` sets each condition's overrides in the
  SBML, resets, and reports which model quantities differ between conditions (`emit.varying`); `write_condition_table_csv`
  writes the independent ones.
* `Experiment.py` loops over that table, one replicate per condition; `Update_parameters.py` gives a replicate its values
  (and sets the initial species directly: RoadRunner evaluates initial assignments only on load);
  `Update_opt_parameters.py` maps a fitted `R_Stage_I_Wuhan` onto the model's `R` in that replicate only.
* PEtab estimates one `sd` per dataset (7); PyAntiGen's default, one profiled sigma per replicate, is the same model.

```powershell
cd problems\Zhao_QuantBiol2020
python Model_run.py --simulate Zhao              # figure in results/Zhao_QuantBiol2020/
python Model_run.py --optimize Zhao --no-fit     # the objective at PEtab's nominal parameters
python Model_run.py --optimize ZhaoAMICI         # WSL: L-BFGS-B with the AMICI gradient
```

Checked against PEtab (`tests/test_zhao.py`): the observables agree to ~4e-9 relative at the nominal and at random parameters;
the benchmark's own `simulations.tsv` is a real simulation here (unlike Zheng's) and the reference reproduces it to 1e-7;
the engine's objective at the nominal parameters is the sum over the seven datasets of `(n/2) ln(SSE/n)` computed from PEtab's
simulation (520.627), with "21 parameters + 7 estimated sigmas = 28", PEtab's own count.

`benchmarks/Zhao_QuantBiol2020/COMPARISON.md` and `COMPRESSION.md` hold the route comparison and the size comparison. The
landscape has several optima, so the three single-start fits end at different points (PEtab's published parameters score NLL
520.6; PyAntiGen on AMICI reaches 387.5, PEtab straight through AMICI 390.7, Nelder-Mead on RoadRunner 426.5): wall-clock
alone does not say which route is better.

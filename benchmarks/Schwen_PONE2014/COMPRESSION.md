# How much does PyAntiGen compress Schwen_PONE2014?

Two questions, measured separately. Sizes: bytes, non-blank lines, **code lines** (no comments or docstrings; for data files, the same as lines) and **gzip -9** bytes, the baseline any representation has to beat to claim it removed information and not only repetition.

## A. What is simulated

The model is small: 14 reactions in PEtab (12 here: the two `Uptake` running totals nobody observes are dropped), two receptor types that follow one pattern (binding, dissociation, internalisation, recycling, degradation), unspecific binding, and observables that are log10 of a fitted scale times a signal. PyAntiGen has one loop over the receptor type and a replicate per condition that sets the insulin the cells start with; the receptor numbers are initial amounts of the fitted `ini_R1` and `ini_R2fold`.

| representation | bytes | lines | code | gzip |
| --- | ---: | ---: | ---: | ---: |
| PEtab SBML (`model.xml`) | 15,732 | 391 | 391 | 1,981 |
| PEtab `conditions.tsv` | 616 | 20 | 20 | 164 |
| **PEtab, what is simulated: SBML + conditions** | **16,348 | 411 | 411 | 2,145** |
| Antimony from Tellurium (`sbmlToAntimony`) | 3,334 | 89 | 89 | 1,201 |
| ... without the `X is "name"` lines, annotations and notes | 2,427 | 58 | 58 | 841 |
| Antimony PyAntiGen generates (12 reactions, rules, parameter and initial values) | 5,447 | 67 | 67 | 1,015 |
| PyAntiGen module: the reactions and rules (`insulin_receptor_Schwen_PONE2014.py`) | 3,987 | 51 | 31 | 829 |
| + value files (`_parameters.csv`, `_InitialConditions.csv`) | 1,298 | 40 | 40 | 640 |
| + conditions table (`data/*_conditions.csv`, 19 rows) | 343 | 20 | 20 | 124 |
| + hooks (`Update_parameters.py`, `Update_opt_parameters.py`) | 1,080 | 14 | 7 | 254 |
| + `Experiment._build_experiment` (the loop over conditions) | 721 | 15 | 15 | 343 |
| **PyAntiGen, what is simulated: all of the above** | **7,429 | 140 | 113 | 2,190** |

In bytes PyAntiGen's definition is 2.2x smaller than PEtab's (7,429 against 16,348); after gzip 1.0x (2,190 against 2,145). Against the Antimony Tellurium writes from the SBML (2,427 bytes clean) the whole PyAntiGen definition is 3.1x as large (113 against 58 code lines); the Antimony it generates is 2.2x as large as the clean one, because the generator writes every species in a rate as `(X / V)` and multiplies by `V` again.

## B. The estimation problem: PEtab tables against Loss_config, Optimizer_settings and the rest

Everything PEtab states besides the SBML and the conditions (which part A counts). `Experiment.py` here also holds the framework's dataclass; the PyAntiGen side contains things PEtab has no place for (the optimizer choice, an AMICI variant, solver tolerances).

| PEtab | bytes | lines | gzip | PyAntiGen | bytes | code lines | gzip | bytes ratio |
| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: |
| parameters (30 estimated: 28 + the two noise parameters; bounds, scale, nominal values, six priors; 30 rows) | 2,183 | 31 | 783 | `Optimizer_settings.py` (nominal values are in the parameters file; the noise parameters are profiled) | 3,996 | 54 | 1,301 | 1.83 |
| observables and noise model (4 observables, log10 scale, 2 noise parameters; the formulas are rules of the model, counted in A) | 942 | 5 | 259 | `Loss_config.py` (which noise parameter each observable shares) | 1,646 | 22 | 451 | 1.75 |
| measurements | 36,095 | 287 | 2,709 | `data/*.csv` + `Data.py` | 11,038 | 282 | 2,952 | 0.31 |
| problem definition and solver | 253 | 13 | 142 | `Experiment.py` + `Solver_settings.py` + two stubs | 4,186 | 67 | 1,455 | 16.55 |
| visualizations | 3,159 | 23 | 403 | `Plots.py` | 2,775 | 49 | 1,044 | 0.88 |
| **total** | **42,632** | **359** | **4,296** | | **23,641** | **474** | **7,203** | **0.55** |

A fixed harness (`AntiGen_paths.py`, `Model_run.py`, `Model_generate.py`: 7,940 bytes, 134 code lines), the same for every model, has no PEtab counterpart.

### The statements themselves

| statement | where | code lines |
| --- | --- | ---: |
| the observables and which noise parameter they share | `Loss_config.py`: `Schwen_loss_config` | 9 |
| the 28 fitted parameters, their bounds, start and scale | `Optimizer_settings.py`: `_build_Schwen_fit` | 17 |
| what distinguishes the 19 conditions | `Update_parameters.py`: `update_Schwen` | 2 |
| one replicate per condition | `Experiment.py`: `_build_experiment` | 15 |
| which data belong to a replicate | `Data.py`: `load_Schwen_data` | 4 |
| the simulated minutes and the output grid | `Solver_settings.py`: `solver_settings_Schwen` | 10 |

Together 57 code lines, plus the reaction module (31 code lines) and two small data files (20 + 272 lines).

## Summary

* **What is simulated:** PyAntiGen's definition is 2.2x smaller in bytes than PEtab's (7,429 against 16,348; 1.0x after gzip). PEtab's SBML is 3.9x the reaction module: 14 kinetic laws as MathML trees, against a loop over the two receptor types. With 12 reactions the model is too small for a loop to save much; the hooks and the conditions table are a large part of PyAntiGen's side.
* **Estimation problem:** 0.55x the bytes of PEtab's tables (23,641 against 42,632). The measurement table is the one place PyAntiGen is much smaller (36,095 bytes against 11,038 for the data file and its loader): PEtab repeats the observable, both condition ids and four observable parameters on every row. The rest is code that PEtab has no place for (the optimizer variants, the loss configuration, the priors as engine settings, the solver grid).
* **Everything, harness excluded:** 31,070 bytes against 58,980 (1.9x).

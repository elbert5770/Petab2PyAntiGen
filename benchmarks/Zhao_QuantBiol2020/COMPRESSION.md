# How much does PyAntiGen compress Zhao_QuantBiol2020?

Two questions, measured separately. Sizes: bytes, non-blank lines, **code lines** (no comments or docstrings; for data files, the same as lines) and **gzip -9** bytes, the baseline any representation has to beat to claim it removed information and not only repetition.

## A. What is simulated

The model is three reactions. PEtab spreads *which region and stage* a simulation is over several places: the six `Trigger_*` parameters of `conditions.tsv`, rules `R := Trigger_Wuhan * (Trigger_Stage_I * R_Stage_I_Wuhan + ...)` (also for gamma_1, gamma_2, the population and alpha, beta), and four initial assignments of the same shape. PyAntiGen has one replicate per (region, stage) that sets its own values.

| representation | bytes | lines | code | gzip |
| --- | ---: | ---: | ---: | ---: |
| PEtab SBML (`model.xml`) | 69,048 | 1,579 | 1,579 | 5,420 |
| PEtab `conditions.tsv` | 507 | 8 | 8 | 152 |
| **PEtab, what is simulated: SBML + conditions** | **69,555 | 1,587 | 1,587 | 5,572** |
| Antimony from Tellurium (`sbmlToAntimony`) | 10,551 | 169 | 169 | 2,869 |
| ... without the `X is "name"` lines, annotations and notes | 7,422 | 123 | 123 | 1,484 |
| Antimony PyAntiGen generates (3 reactions, rules, parameter and initial values) | 1,663 | 21 | 21 | 479 |
| PyAntiGen module: the reactions and rules (`siqc_Zhao_QuantBiol2020.py`) | 1,823 | 24 | 13 | 422 |
| + value files (`_parameters.csv`, `_InitialConditions.csv`) | 678 | 12 | 12 | 418 |
| + conditions table (`data/*_conditions.csv`, 7 rows) | 599 | 8 | 8 | 268 |
| + hooks (`Update_parameters.py`, `Update_opt_parameters.py`) | 1,681 | 25 | 14 | 412 |
| + `Experiment._build_experiment` (the loop over conditions) | 771 | 15 | 15 | 383 |
| **PyAntiGen, what is simulated: all of the above** | **5,552 | 84 | 62 | 1,903** |

In bytes PyAntiGen's definition is 12.5x smaller than PEtab's (5,552 against 69,555); after gzip 2.9x (1,903 against 5,572). Against the Antimony Tellurium writes from the SBML (7,422 bytes clean), the whole PyAntiGen definition is 1.3x smaller, and that Antimony still carries all 37 parameters and the trigger rules.

## B. The estimation problem: PEtab tables against Loss_config, Optimizer_settings and the rest

Everything PEtab states besides the SBML and the conditions (which part A counts). `Experiment.py` here also holds the framework's dataclass; the PyAntiGen side contains things PEtab has no place for (the optimizer choice, an AMICI variant, solver tolerances).

| PEtab | bytes | lines | gzip | PyAntiGen | bytes | code lines | gzip | bytes ratio |
| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: |
| parameters (28 estimated: 21 rates, 7 noise; bounds, scale, nominal values) | 2,380 | 30 | 458 | `Optimizer_settings.py` (nominal values are in the conditions table) | 2,477 | 39 | 769 | 1.04 |
| observables and noise model (5 observables, 1 measured) | 682 | 6 | 191 | `Loss_config.py` (one observable; one sigma per replicate is the default) | 591 | 9 | 189 | 0.87 |
| measurements | 12,813 | 83 | 722 | `data/*.csv` + `Data.py` | 3,118 | 94 | 763 | 0.24 |
| problem definition and solver | 253 | 13 | 142 | `Experiment.py` + `Solver_settings.py` + two stubs | 3,445 | 57 | 1,228 | 13.62 |
| visualizations | 1,700 | 8 | 319 | `Plots.py` | 1,214 | 23 | 557 | 0.71 |
| **total** | **17,828** | **140** | **1,832** | | **10,845** | **222** | **3,506** | **0.61** |

A fixed harness (`AntiGen_paths.py`, `Model_run.py`, `Model_generate.py`: 7,758 bytes, 132 code lines), the same for every model, has no PEtab counterpart.

### The statements themselves

| statement | where | code lines |
| --- | --- | ---: |
| the observable and its noise | `Loss_config.py`: `Zhao_loss_config` | 7 |
| the 21 fitted parameters, their bounds, start and scale | `Optimizer_settings.py`: `_build_Zhao_fit` | 19 |
| what distinguishes the 7 conditions, and how a fitted name reaches the model | `Update_parameters.py`: `update_Zhao`; `Update_opt_parameters.py`: `update_opt_Zhao` | 10 |
| one replicate per condition | `Experiment.py`: `_build_experiment` | 15 |
| which data belong to a replicate | `Data.py`: `load_Zhao_data` | 3 |
| the simulated days | `Solver_settings.py`: `solver_settings_Zhao` | 2 |

Together 56 code lines, plus the reaction module (13 code lines) and two small data files (8 + 83 lines).

## Summary

* **What is simulated:** 12.5x smaller in bytes (2.9x after gzip). PEtab's SBML is 38x the reaction module because it also carries the trigger rules and initial assignments that PyAntiGen replaces with a seven-row table and two hooks.
* **Estimation problem:** 0.61x the bytes of PEtab's tables (10,845 against 17,828), and almost all of the saving is the measurement table: PEtab repeats the observable, the condition and the noise parameter on each of its 82 rows (12,813 bytes against 3,118 for the data file and its loader). The parameters and the problem definition are larger in PyAntiGen, which carries framework structure there.
* **Everything, harness excluded:** 16,397 bytes against 87,383 (5.3x).

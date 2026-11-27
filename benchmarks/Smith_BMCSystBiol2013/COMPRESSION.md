# How much does PyAntiGen compress Smith_BMCSystBiol2013?

Two questions, measured separately. Sizes: bytes, non-blank lines, **code lines** (no comments or docstrings; for data files, the same as lines) and **gzip -9** bytes, the baseline any representation has to beat to claim it removed information and not only repetition.

## A. What is simulated

The model is 313 reactions after merging, most of them the FOXO1 network: 16 modification states (Pa, Pd, Pe, pUb) in cytoplasm, nucleus and DNA-bound, with transport, phosphorylation, ubiquitination and degradation. PEtab writes each of its 367 reactions as an SBML kinetic law with its own rate parameter (about 150 `R<n>_ktr`, `R<n>_kkin`, ... values that are twelve distinct numbers). PyAntiGen has the loop and a base rate and a factor per modification, and one replicate per condition that sets Ins, H2O2, the two on/off switches and E2F1.

| representation | bytes | lines | code | gzip |
| --- | ---: | ---: | ---: | ---: |
| PEtab SBML (`model.xml`) | 350,061 | 8,871 | 8,871 | 13,618 |
| PEtab `conditions.tsv` | 2,694 | 36 | 36 | 441 |
| **PEtab, what is simulated: SBML + conditions** | **352,755 | 8,907 | 8,907 | 14,059** |
| Antimony from Tellurium (`sbmlToAntimony`) | 105,504 | 1,347 | 1,347 | 12,435 |
| ... without the `X is "name"` lines, annotations and notes | 76,905 | 1,037 | 1,037 | 10,114 |
| Antimony PyAntiGen generates (313 reactions, rules, parameter and initial values) | 73,125 | 628 | 628 | 6,117 |
| PyAntiGen module: the reactions and rules (`insulin_ros_foxo_Smith_BMCSystBiol2013.py`) | 15,961 | 192 | 152 | 3,237 |
| + value files (`_parameters.csv`, `_InitialConditions.csv`) | 5,218 | 212 | 212 | 1,447 |
| + conditions table (`data/*_conditions.csv`, 36 rows) | 1,413 | 36 | 36 | 356 |
| + hook (`Update_parameters.py`) | 624 | 8 | 4 | 158 |
| + `Experiment._build_experiment` (the loop over conditions) | 725 | 15 | 15 | 345 |
| **PyAntiGen, what is simulated: all of the above** | **23,941 | 463 | 419 | 5,543** |

In bytes PyAntiGen's definition is 14.7x smaller than PEtab's (23,941 against 352,755); after gzip 2.5x (5,543 against 14,059). Against the Antimony Tellurium writes from the SBML (76,905 bytes clean) the whole PyAntiGen definition is 3.2x smaller (419 against 1037 code lines); the Antimony it generates is 1.0x as large as the clean one, because the generator writes every species as `V * (X / V)` here (the laws are in amounts and PyAntiGen reads a species as a concentration) and the 48 FOXO1 species are spelled out.

## B. The estimation problem: PEtab tables against Loss_config, Optimizer_settings and the rest

Everything PEtab states besides the SBML and the conditions (which part A counts). `Experiment.py` here also holds the framework's dataclass; the PyAntiGen side contains things PEtab has no place for (the optimizer choice, an AMICI variant, solver tolerances).

| PEtab | bytes | lines | gzip | PyAntiGen | bytes | code lines | gzip | bytes ratio |
| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: |
| parameters (25 estimated: 16 rate constants, 9 observable scales; bounds, scale, nominal values; 80 rows, 55 fixed) | 3,402 | 80 | 1,094 | `Optimizer_settings.py` (nominal values are in the parameters file) | 4,078 | 48 | 1,063 | 1.20 |
| observables and noise model (9 observables, sigma 1; the formulas are rules of the model, counted in A) | 494 | 10 | 273 | `Loss_config.py` (the declared sigma) | 1,174 | 21 | 405 | 2.38 |
| measurements | 3,532 | 63 | 642 | `data/*.csv` + `Data.py` | 2,758 | 61 | 785 | 0.78 |
| problem definition and solver | 253 | 13 | 142 | `Experiment.py` + `Solver_settings.py` + two stubs | 4,399 | 65 | 1,440 | 17.39 |
| visualizations | 3,317 | 24 | 553 | `Plots.py` | 2,999 | 48 | 1,163 | 0.90 |
| **total** | **10,998** | **190** | **2,704** | | **15,408** | **243** | **4,856** | **1.40** |

A fixed harness (`AntiGen_paths.py`, `Model_run.py`, `Model_generate.py`: 8,613 bytes, 145 code lines), the same for every model, has no PEtab counterpart.

### The statements themselves

| statement | where | code lines |
| --- | --- | ---: |
| the observables and their (declared) noise | `Loss_config.py`: `Smith_loss_config` | 10 |
| the 25 fitted parameters, their bounds, start and scale | `Optimizer_settings.py`: `_build_Smith_fit` | 17 |
| what distinguishes the 36 conditions | `Update_parameters.py`: `update_Smith` | 3 |
| one replicate per condition | `Experiment.py`: `_build_experiment` | 15 |
| which data belong to a replicate | `Data.py`: `load_Smith_data` | 4 |
| the simulated minutes and the output grid | `Solver_settings.py`: `solver_settings_Smith` | 8 |

Together 57 code lines, plus the reaction module (152 code lines) and two small data files (36 + 51 lines).

## Summary

* **What is simulated:** 14.7x smaller in bytes (2.5x after gzip). PEtab's SBML is 22x the reaction module: 367 kinetic laws, each with its MathML tree and (for FOXO1) its own rate parameter, against loops over the modification states.
* **Estimation problem:** 1.40x the bytes of PEtab's tables (15,408 against 10,998): PyAntiGen is smaller only where PEtab repeats itself, the measurement table (3,532 bytes against 2,758 for the data file and its loader) and the plots. It is larger for the parameters (the three optimizer variants and the solver notes; the 55 fixed parameters are in the value file, counted in A), the noise model (one declaration per observable in code against ten table rows) and for the problem definition (the dataclass, the loop over conditions and the output grid, which PEtab leaves to the simulator). Here the model is the saving, not the estimation problem.
* **Everything, harness excluded:** 39,349 bytes against 363,753 (9.2x).

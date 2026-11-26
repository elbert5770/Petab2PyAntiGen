# How much does PyAntiGen compress Sneyd_PNAS2002?

Two questions, measured separately. Sizes: bytes, non-blank lines, **code lines** (no comments or docstrings; for data files, the same as lines) and **gzip -9** bytes, the baseline any representation has to beat to claim it removed information and not only repetition.

## A. What is simulated

The model is ten first-order transitions between six receptor states, whose rate coefficients depend on the calcium and IP3 of the experiment. PEtab states each condition as a row of `conditions.tsv` (two numbers) and writes each reaction as a full SBML kinetic law; PyAntiGen has one replicate per condition that sets its IP3 and Ca2+, and the ten transitions as a list of (from, to, coefficient) that one loop turns into reactions.

| representation | bytes | lines | code | gzip |
| --- | ---: | ---: | ---: | ---: |
| PEtab SBML (`model.xml`) | 16,993 | 483 | 483 | 2,177 |
| PEtab `conditions.tsv` | 595 | 10 | 10 | 178 |
| **PEtab, what is simulated: SBML + conditions** | **17,588 | 493 | 493 | 2,355** |
| Antimony from Tellurium (`sbmlToAntimony`) | 2,529 | 86 | 86 | 1,110 |
| ... without the `X is "name"` lines, annotations and notes | 1,987 | 67 | 67 | 842 |
| Antimony PyAntiGen generates (10 reactions, rules, parameter and initial values) | 2,436 | 44 | 44 | 728 |
| PyAntiGen module: the transitions and rules (`ip3r_Sneyd_PNAS2002.py`) | 2,444 | 41 | 27 | 606 |
| + value files (`_parameters.csv`, `_InitialConditions.csv`) | 684 | 25 | 25 | 426 |
| + conditions table (`data/*_conditions.csv`, 9 rows) | 256 | 10 | 10 | 109 |
| + hook (`Update_parameters.py`) | 450 | 6 | 3 | 76 |
| + `Experiment._build_experiment` (the loop over conditions) | 739 | 15 | 15 | 354 |
| **PyAntiGen, what is simulated: all of the above** | **4,573 | 97 | 80 | 1,571** |

In bytes PyAntiGen's definition is 3.8x smaller than PEtab's (4,573 against 17,588); after gzip 1.5x (1,571 against 2,355). Against the Antimony Tellurium writes from the SBML (1,987 bytes clean) the whole PyAntiGen definition is 2.3x as large (80 against 67 code lines), and the Antimony it generates is 1.2x as large as the clean one: the generated Antimony divides by the explicit volume `V_Membrane` and multiplies by it again in every rate. What PyAntiGen saves over PEtab here is the SBML's verbosity (ids, annotations, units, MathML), most of which Tellurium also drops.

## B. The estimation problem: PEtab tables against Loss_config, Optimizer_settings and the rest

Everything PEtab states besides the SBML and the conditions (which part A counts). `Experiment.py` here also holds the framework's dataclass; the PyAntiGen side contains things PEtab has no place for (the optimizer choice, an AMICI variant, solver tolerances).

| PEtab | bytes | lines | gzip | PyAntiGen | bytes | code lines | gzip | bytes ratio |
| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: |
| parameters (15 estimated: 14 rate constants, 1 noise; bounds, scale, nominal values) | 772 | 16 | 333 | `Optimizer_settings.py` (nominal values are in the parameters file) | 2,879 | 45 | 920 | 3.73 |
| observables and noise model (1 observable; its formula is a rule of the model, counted in A) | 194 | 2 | 148 | `Loss_config.py` (one observable; the shared sigma is the `sigma_block` in `Optimizer_settings.py`) | 1,218 | 15 | 290 | 6.28 |
| measurements | 15,863 | 136 | 934 | `data/*.csv` + `Data.py` | 4,800 | 145 | 898 | 0.30 |
| problem definition and solver | 253 | 13 | 142 | `Experiment.py` + `Solver_settings.py` + two stubs | 3,494 | 57 | 1,207 | 13.81 |
| visualizations | 1,557 | 10 | 288 | `Plots.py` | 1,243 | 24 | 578 | 0.80 |
| **total** | **18,639** | **177** | **1,845** | | **13,634** | **286** | **3,893** | **0.73** |

A fixed harness (`AntiGen_paths.py`, `Model_run.py`, `Model_generate.py`: 8,134 bytes, 139 code lines), the same for every model, has no PEtab counterpart.

### The statements themselves

| statement | where | code lines |
| --- | --- | ---: |
| the observable and its noise | `Loss_config.py`: `Sneyd_loss_config` | 7 |
| the 14 fitted rate constants, their bounds, start and scale (one shared sigma) | `Optimizer_settings.py`: `_build_Sneyd_fit` | 16 |
| what distinguishes the 9 conditions | `Update_parameters.py`: `update_Sneyd` | 3 |
| one replicate per condition | `Experiment.py`: `_build_experiment` | 15 |
| which data belong to a replicate | `Data.py`: `load_Sneyd_data` | 3 |
| the simulated seconds | `Solver_settings.py`: `solver_settings_Sneyd` | 2 |

Together 46 code lines, plus the reaction module (27 code lines) and two small data files (10 + 136 lines).

## Summary

* **What is simulated:** 3.8x smaller in bytes (1.5x after gzip). PEtab's SBML is 7.0x the reaction module: ten kinetic laws, each with its MathML tree, against ten one-line coefficients.
* **Estimation problem:** 0.73x the bytes of PEtab's tables (13,634 against 18,639), and almost all of the saving is the measurement table: PEtab repeats the observable, the condition, the noise parameter and four bookkeeping columns on each of its 135 rows (15,863 bytes against 4,800 for the data file and its loader). The parameters and the problem definition are larger in PyAntiGen, which carries framework structure there.
* **Everything, harness excluded:** 18,207 bytes against 36,227 (2.0x).

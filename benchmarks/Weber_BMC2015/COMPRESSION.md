# How much does PyAntiGen compress Weber_BMC2015?

Two questions, measured separately. Sizes: bytes, non-blank lines, **code lines** (no comments or docstrings; for data files, the same as lines) and **gzip -9** bytes, the baseline any representation has to beat to claim it removed information and not only repetition.

## A. What is simulated

The model is 20 reactions. PEtab says what a condition does to it with six condition parameters per condition (an on/off flag, a dose and a time for each stimulus), which SBML rules such as `u5 := piecewise(0, time - PdBu_time < 0, PdBu_dose)` turn into the four inputs. PyAntiGen has the inputs as parameters that an event sets, so the whole of `conditions.tsv` is four lines in `Events.py`. The observables, which PEtab keeps in `observables.tsv`, are rules of the module (a signal times a fitted scale factor).

| representation | bytes | lines | code | gzip |
| --- | ---: | ---: | ---: | ---: |
| PEtab SBML (`model.xml`) | 22,635 | 623 | 623 | 2,610 |
| PEtab `conditions.tsv` | 204 | 4 | 4 | 140 |
| **PEtab, what is simulated: SBML + conditions** | **22,839 | 627 | 627 | 2,750** |
| Antimony from Tellurium (`sbmlToAntimony`) | 4,064 | 126 | 126 | 1,648 |
| ... without the `X is "name"` lines, annotations and notes | 3,038 | 86 | 86 | 1,237 |
| Antimony PyAntiGen generates (20 reactions, rules, parameter and initial values) | 3,877 | 79 | 79 | 1,182 |
| PyAntiGen module: the reactions, the observable rules (`pkd_pi4k_cert_Weber_BMC2015.py`) | 4,143 | 55 | 36 | 841 |
| + value files (`_parameters.csv`, `_InitialConditions.csv`) | 1,392 | 45 | 45 | 712 |
| + the conditions, as events (`Events.py`) | 879 | 13 | 10 | 338 |
| + `Experiment._build_experiment` (one replicate per condition) | 590 | 13 | 13 | 266 |
| **PyAntiGen, what is simulated: all of the above** | **7,004 | 126 | 104 | 2,157** |

In bytes PyAntiGen's definition is 3.3x smaller than PEtab's (7,004 against 22,839); after gzip 1.3x (2,157 against 2,750). Against the Antimony Tellurium writes from the SBML (3,038 bytes clean) the whole PyAntiGen definition is 2.3x as large (104 against 86 code lines), and the Antimony it generates is 1.3x as large as the clean one. There is nothing to loop over (20 distinct reactions), the generated Antimony divides by the explicit volume `V_Cyt` and multiplies by it again in every rate, and the definition carries the experiment scaffolding. What PyAntiGen saves over PEtab here is the SBML's verbosity (ids, annotations, MathML, the condition rules), most of which Tellurium also drops, and not the looping that made Zheng compress.

## B. The estimation problem: PEtab tables against Loss_config, Optimizer_settings and the rest

Everything PEtab states besides the SBML and the conditions (which part A counts). `Experiment.py` here also holds the framework's dataclass; the PyAntiGen side contains things PEtab has no place for (the optimizer choice, an AMICI variant, solver tolerances).

| PEtab | bytes | lines | gzip | PyAntiGen | bytes | code lines | gzip | bytes ratio |
| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: |
| parameters (36 estimated: 31 rates, inputs and scales, 5 noise; 4 fixed; bounds, scale, nominal values) | 2,136 | 41 | 772 | `Optimizer_settings.py` (nominal values are in the parameters file) | 2,993 | 48 | 1,031 | 1.40 |
| observables and noise model (8 observables; formulas are in the module, counted in A) | 780 | 9 | 248 | `Loss_config.py` (the three declared sigmas; the others are profiled by default) | 1,517 | 20 | 482 | 1.94 |
| measurements | 11,169 | 136 | 1,504 | `data/*.csv` + `Data.py` | 4,014 | 83 | 1,308 | 0.36 |
| problem definition and solver | 253 | 13 | 142 | `Experiment.py` + `Solver_settings.py` + `Observed_species.py` (+ `Events.py`, counted in A) | 3,767 | 55 | 1,063 | 14.89 |
| visualizations | 710 | 9 | 221 | `Plots.py` | 1,654 | 31 | 718 | 2.33 |
| **total** | **15,048** | **208** | **2,887** | | **13,945** | **237** | **4,602** | **0.93** |

A fixed harness (`AntiGen_paths.py`, `Model_run.py`, `Model_generate.py`: 7,872 bytes, 134 code lines), the same for every model, has no PEtab counterpart.

### The statements themselves

| statement | where | code lines |
| --- | --- | ---: |
| the observables and their noise (declared or profiled) | `Loss_config.py`: `Weber_loss_config`, `measured` | 13 |
| the 31 fitted parameters, their bounds, start and scale | `Optimizer_settings.py`: `_build_Weber_fit` | 15 |
| what each condition does (the inputs, and when) | `Events.py`: `Weber_events` | 3 |
| one replicate per condition | `Experiment.py`: `_build_experiment` | 13 |
| which data belong to a replicate | `Data.py`: `load_Weber_data` | 5 |
| pre-equilibration and the observed hours | `Solver_settings.py`: `solver_settings_Weber` | 5 |

Together 54 code lines, plus the reaction module (36 code lines) and one data file (75 lines).

## Summary

* **What is simulated:** 3.3x smaller in bytes (1.3x after gzip). PEtab's SBML is 5.5x the reaction module; the module also carries the observables, which PEtab keeps in a separate table.
* **Estimation problem:** 0.93x the bytes of PEtab's tables (13,945 against 15,048), and almost all of the saving is the measurement table: PEtab repeats the observable, both conditions and the scale and noise parameters on each of its 135 rows (11,169 bytes against 4,014 for the data file and its loader, where replicates share a row). The parameters and the problem definition are larger in PyAntiGen, which carries framework structure there.
* **Everything, harness excluded:** 20,949 bytes against 37,887 (1.8x).

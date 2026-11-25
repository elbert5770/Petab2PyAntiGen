# How much does PyAntiGen compress Zheng_PNAS2012?

Two questions, measured separately. Sizes: bytes, non-blank lines, **code lines** (no comments or docstrings; for data files, the same as lines) and **gzip -9** bytes, the baseline any representation has to beat to claim it removed information and not only repetition.

## A. The model: how much text defines the reactions?

| representation | bytes | lines | code | gzip |
| --- | ---: | ---: | ---: | ---: |
| PEtab SBML (`model.xml`, one function definition per reaction) | 92,169 | 2,832 | 2,832 | 4,733 |
| SBML exported from PyAntiGen's Antimony (for scale) | 51,945 | 1,361 | 1,361 | 1,857 |
| Antimony from Tellurium (`sbmlToAntimony`) | 9,617 | 221 | 221 | 2,512 |
| ... without the `X is "name"` lines, annotations and notes | 7,877 | 148 | 148 | 1,876 |
| Antimony PyAntiGen generates (reactions, rules, parameter and initial values) | 10,322 | 139 | 139 | 1,098 |
| **PyAntiGen module: lattice, reaction loops** (`methylation_Zheng_PNAS2012.py`) | 2,953 | 49 | 30 | 592 |
| module + the value files it needs (`_parameters.csv`, `_InitialConditions.csv`) | 5,046 | 113 | 94 | 1,537 |

Compression of the *definition* (rows against the module alone, ratio of that representation to the module):

| from | bytes | code lines | gzip |
| --- | ---: | ---: | ---: |
| PEtab SBML | 31.2x | n/a | 8.0x |
| Tellurium Antimony | 3.3x | 7.4x | 4.2x |
| clean Antimony | 2.7x | 4.9x | 3.2x |
| PyAntiGen's generated Antimony | 3.5x | 4.6x | 1.9x |

Reading it honestly: the Antimony PyAntiGen generates is **not** shorter than Tellurium's, it is the same network written out (`k * (S/V) * V` with the volume explicit). The gain is in what a person has to write and read: the generator, not its output. gzip shows how much of the SBML is mere repetition: it shrinks 19x on its own, so the module's 8.0x advantage over the gzipped SBML is the part that is more than repetition. Against the generated Antimony the gzip ratio is 1.9x: that file is itself very regular, and the module is little more than its irreducible content. (Bytes are counted with LF line endings, so they do not depend on the platform; files on a Windows checkout are a little larger.)

### Growth with the size of the system

The same module, `N_ME` changed from 4 to other lattice sizes (an n x n lattice of methylation states): its text does not change; the Antimony and SBML it stands for grow with the number of reactions.

| lattice | states | reactions | module code lines | generated Antimony (bytes) | SBML from it (bytes) | SBML / module |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2 x 2 | 3 | 8 | 30 | 1,423 | 7,578 | 3x |
| 3 x 3 | 8 | 29 | 30 | 4,992 | 25,456 | 9x |
| 4 x 4 | 15 | 60 | 30 | 10,322 | 51,945 | 18x |
| 6 x 6 | 35 | 152 | 30 | 26,018 | 130,509 | 44x |
| 8 x 8 | 63 | 284 | 30 | 48,451 | 243,210 | 82x |
| 12 x 12 | 143 | 668 | 30 | 115,621 | 572,978 | 194x |

(Reactions counted as the engine writes them: each reversible step is two.)

## B. The estimation problem: PEtab tables against Experiment, Loss_config and the rest

Everything PEtab states besides the SBML, against the PyAntiGen files that say the same. The PyAntiGen side also contains things PEtab has no place for (three optimizer specs, solver and backend settings, a known-sigma loss variant for Levenberg-Marquardt) and the observables' *formulas*, which PEtab keeps in `observables.tsv`.

| PEtab | bytes | lines | gzip | PyAntiGen | bytes | code lines | gzip | bytes ratio |
| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: |
| parameters (nominal values, bounds, scale, estimate) | 2,659 | 47 | 791 | `Optimizer_settings.py` + `_parameters.csv` | 4,324 | 98 | 1,543 | 1.63 |
| observables and noise model | 5,187 | 16 | 327 | `Loss_config.py` + the `Observables` class | 1,943 | 27 | 649 | 0.37 |
| conditions, pre-equilibration | 71 | 3 | 74 | `Events.py` + `Solver_settings.py` | 2,052 | 25 | 567 | 28.90 |
| measurements | 4,587 | 61 | 523 | `data/*.csv` + `Data.py` | 1,457 | 14 | 489 | 0.32 |
| problem definition (what is simulated, what is fitted) | 253 | 13 | 142 | `Experiment.py` + three small hooks | 2,621 | 47 | 846 | 10.36 |
| visualizations | 2,037 | 16 | 328 | `Plots.py` | 1,326 | 27 | 633 | 0.65 |
| **total** | **14,794** | **156** | **2,185** | | **13,723** | **238** | **4,727** | **0.93** |

On top of that, PyAntiGen needs a fixed harness (`AntiGen_paths.py`, `Model_run.py`, `Model_generate.py`: 8,130 bytes, 140 code lines) that is the same for every model and has no PEtab counterpart.

### What is information and what is boilerplate

A good deal of PyAntiGen's problem statement is structure that a model does not choose: the `Optimization` and `Experiment` dataclasses, the hooks that do nothing here (`Update_parameters`, `Update_opt_parameters`, `Observed_species`), and the plot code. The lines that *state the problem*:

| statement | where | code lines |
| --- | --- | ---: |
| the 15 observables and their noise model | `Loss_config.py`: `Zheng_loss_config` | 9 |
| the observables' formulas (fractions of total H3) | the `Observables` class | 7 |
| which parameters are fitted, from where, on what scale, within what bounds | `Optimizer_settings.py`: `_build_Zheng_fit` | 21 |
| the experiment: pre-equilibrate, switch dilution off, observe 25 h | `Solver_settings.py`: `solver_settings_Zheng`; `Events.py`: `Zheng_event` | 7 |
| which data belong to the experiment | `Data.py`: `load_Zheng_data` | 4 |

Together 48 code lines, against 156 lines in PEtab's tables. The rest of PyAntiGen's 238 code lines is framework (the `Optimization` and `Experiment` dataclasses, the solver factory, backend settings, plot code, empty hooks) and things PEtab cannot say at all (optimizer choice, solver tolerances, a second loss for Levenberg-Marquardt).

## Summary

* **Model.** The module is 31x smaller than PEtab's SBML in bytes (8x after gzip, which already removes the SBML's repetition) and 3.5x smaller than the Antimony it generates (1.9x after gzip). The generated Antimony is as long as Tellurium's: nothing is gained in the Antimony itself, only in the text a person writes. The module does not change as the system grows while the SBML grows with the reactions (194x the module at a 12 x 12 lattice).
* **Estimation problem.** In bytes it is 0.93x PEtab's tables (13,723 against 14,794): no real saving, because PEtab's tables are compact for one condition and 15 observables and PyAntiGen's problem files carry framework structure. The statements themselves are 48 code lines against 156 table lines: the 15 observables are one loop, the fitted parameters come from the model's own rate-constant helper, and a different sigma model or optimizer is a one-key change. Conditions are where PEtab is clearly more compact (one 71-byte table against the solver blocks and an event).

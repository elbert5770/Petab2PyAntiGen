---
name: Convert a PEtab problem to PyAntiGen
description: Workflow, conventions and pitfalls for replacing a PEtab v1 problem (SBML + tables) with a PyAntiGen project, verified against PEtab. Worked example problems/Zheng_PNAS2012.
---

# Convert a PEtab problem to PyAntiGen

Read `problems/Zheng_PNAS2012/` first: it is the worked example, and
`python -m petab2pyantigen stats ...` shows what the conversion saves.

## 1. Look at the model, not at the SBML

`python -m petab2pyantigen antimony <model.xml>` gives clean Antimony (Tellurium, with
`X is "name"` lines, annotations and notes removed; `comp` modules are flattened).
Read it for *structure*: which indices do the species and rate constants share
(`K27me<a>K36me<b>`, `k<ab>_<cd>`)? Which reactions are one pattern repeated? A
SBML function definition per reaction (`Function_for_v_17`) is PEtab-exporter noise; the
rate law is just the body. Do not translate reaction by reaction.

## 2. Conventions (PyAntiGen)

* species = `<species>_<Compartment>`; with isotopes `<species>_<isotope>_<Compartment>`.
* compartment names CamelCase, **no underscores** (SBML `default` -> `Default`);
  every compartment has an explicit `V_<Compartment>` parameter.
* species are substance-only (amounts); a rate law `k * [S] * V` is `MA` with prototype `k`.
* `RMA` with `"[k_forward, k_reverse]"` writes both directions of a reversible step.
* empty reactants / products: pass `""`.
* never hard-code numbers in a module; use parameter names. Values go in
  `antimony_models/<Model>/<Model>_parameters.csv`.
* keep PEtab's parameter and observable ids (`k00_01`, `observable_K27me0K36me0`), so every
  PEtab table row has an obvious PyAntiGen counterpart.
* one module file + class per problem, named after it (`methylation_<Model>.py`).

## 3. Where each PEtab piece goes

| PEtab | PyAntiGen |
| --- | --- |
| SBML reactions, rules | `antimony_modules/<Area>/<module>.py` (f-string loops, `add_reaction`, `add_rule`) |
| parameters.tsv nominalValue | `<Model>_parameters.csv` (`python -m petab2pyantigen emit`) |
| parameters.tsv bounds, scale, estimate | `Modules/Optimizer_settings.py` (`param_names`, `x0`, `bounds`, `parameter_scale`) |
| SBML initial values | `<Model>_InitialConditions.csv` (amounts, not concentrations) |
| conditions.tsv | `Modules/Events.py` (switches), `Update_parameters.py` (per-replicate values) |
| preequilibrationConditionId | an untracked first block in `Solver_settings.py`; shift data times by its length in `Data.py` |
| observables.tsv | assignment rules in the module + `Modules/Loss_config.py` |
| noiseFormula / noiseParameters | `"sigma_pool"` / `"sigma_block"` / `sigma_method: "fixed"` in the loss_config: see pitfalls |
| measurements.tsv | `data/<Model>.csv`, wide: one row per time, one column per observable |
| visualizations.tsv | `Modules/Plots.py` |

## 4. Verify -- always

`python -m petab2pyantigen verify --petab <dir> --project problems/<Model> --experiment <key>`
simulates both sides (PEtab: SBML + RoadRunner steady state; PyAntiGen: its engine) at the
nominal parameters and at random ones, and compares every observable and the log-likelihood.
Expect agreement of ~1e-8 or better. `tests/test_zheng.py` is the template for a regression test.

## 5. Pitfalls found on Zheng_PNAS2012

* **Do not trust `simulations.tsv`.** In the benchmark collection it can simply be a copy of the
  measurements (Zheng's is: chi2 ~ 1e-26). Check that the reference you compare with is real:
  at the nominal parameters chi2 == N if sigma was fitted (the MLE).
* **Pre-equilibration** to steady state is a long untracked block (relaxation time is
  1/(smallest eigenvalue); here dilution makes it ~ 1/inflowp). Switch the condition with an
  event at the block end and add that time to the data times.
* **One shared sigma.** By default the engine estimates one sigma per observable, which for
  15 observables x 4 points warns and over-fits. PEtab's single `sigma` in every `noiseFormula`
  is `"sigma_pool": "sigma"` on every observable entry of the `loss_config`
  (`Loss_config.py`, `SHARED_SIGMA`). Needs a PyAntiGen with `sigma_pool` (see requirements.txt).
  `sigma_block` on the loss element pools one observable across simulations; the two compose.
* **Import order in `Model_run.py`:** `import AntiGen_paths` before `Modules`, because the
  modules import `antimony_modules` and `AntiGen_paths` puts the repository root on `sys.path`.
* **libsbml:** a `Model` dangles once its `SBMLDocument` is garbage collected; keep the document.
* **AMICI (WSL):** put `/opt/venv/bin` first on `PATH` (SWIG 4.5.1 must match the one AMICI was built with;
  the system SWIG 4.2.0 gives "SwigPyObject has no attribute module" at load), and the PyAntiGen checkout on
  `PYTHONPATH`. Nelder-Mead ignores `analytic_gradient`; use L-BFGS-B (LM additionally needs known sigmas).
  A spec with `analytic_gradient` runs on AMICI without `--backend`.
  Levenberg-Marquardt (`ZhengLM`) is least squares, so every sigma must be declared (`sigma_method: "fixed"`); with a
  profiled or pooled sigma it stops at "no evaluable start". L-BFGS-B keeps the concentrated likelihood.
  Observables that are bare result columns (`observable_K27me0K36me0`) get an exact chain-rule gradient; an
  expression (`SF*pred`) or a callable observable falls back to ~600 finite-difference loss evaluations per gradient
  (~1.9 s on Zheng). Prefer assignment rules in the model over expressions in the loss config.
* **Direct AMICI as a yardstick:** `benchmarks/amici_direct.py` runs the untouched PEtab problem through AMICI
  (59 s for the Zheng fit); `benchmarks/<Problem>/run_all_fits.sh` + `compare_fits.py --problem <Problem>` score every route with PEtab's own objective. PEtab's explicit `sigma` must start near its ML value or the joint fit stalls; PyAntiGen's
  profiled (concentrated) sigma has no such problem.
* **`six`** is needed by Tellurium but no longer declared by it (it is in requirements.txt).

## 6. Pitfalls found on Zhao_QuantBiol2020 (conditions that are switches)

* **Conditions as 0/1 triggers.** PEtab can pick each condition's parameters and initial states with indicator parameters
  multiplied into SBML rules. Do not translate the rules: `emit.effective_condition_values(problem)` evaluates what each
  condition does (set its overrides, `reset()`), `emit.varying` lists what differs, and the independent columns become a small
  per-condition CSV that `Experiment.py` loops over, one replicate per condition.
* **Initial states in `Update_parameters` must be set on the species.** RoadRunner evaluates a species' initial assignment
  when the model loads, and the engine's order is reset -> `Update_parameters` -> fitted values -> hooks, so a parameter
  changed there does not reach an initial assignment that reads it.
* **A fitted name only a hook understands** (`R_Stage_I_Wuhan` -> `R` of one replicate): `Update_opt_parameters`. `verify.py`
  applies overrides the way a fit does, so this works there too; PyAntiGen's `member=` refuses a parameter a replicate does
  not read, which is why verify does not use it.
* **AMICI needs `"no_inputs": True`** on a replicate with no events (nothing clock-driven), or it asks for an input schedule.
* **AMICI used to see stale values** for plain parameters `Update_parameters` sets but no fit moves (a regional population),
  and to refuse data at the window's first time (day 0): the fit silently fell back to finite differences ("the analytic
  gradient failed N time(s)": look for that line). Both are fixed in PyAntiGen; check the log of any AMICI fit for it.
* **The reference must apply the first condition before the initial state is evaluated** (`reference.simulate_petab` does).
* **Several optima:** one start per route ends at different points; compare objective values, not only times.

## 7. Pitfalls found on Weber_BMC2015 (scale parameters, fixed noise, clock-driven inputs)

* **Observable scale parameters** (`observableFormula = signal * observableParameter1_y`, with `observableParameters` naming a
  fitted `scale_y`): write the observable as an assignment rule of the module (`yPKDpN0 := PKDDAGa_Cyt * scale_yPKDpN0`) and make
  `scale_y` an ordinary model parameter. It is fitted by name, with no hook. `reference.simulate_petab` substitutes the
  `observableParameter<i>_<obs>` placeholders from each measurement row.
* **Noise parameters that PEtab does not estimate** (`estimate = 0`, here the three "total" observables, with sigmas of 1e6 to 2e8):
  declare them in `Loss_config` with `"sigma_method": "fixed", "sigma_value": ...`. The estimated ones stay profiled. The two
  conventions differ by a constant (profiled `(n/2) ln(SSE/n)`, declared `n ln(sigma) + SSE/(2 sigma^2)`), and
  `benchmarks/compare_fits.py` scores each group the way its noise parameter is treated.
* **One noise parameter for an observable measured in several conditions** (`std_yPKDpN0`): the loss elements share a
  `"sigma_block": "<name>"`, which makes one sigma per observable across them. `sigma_pool` pools observables inside one simulation
  and does not do this: two blocks appear, the objective is off by ~10, and nothing fails. Check the engine's line
  "N block(s), M points, k parameters + S estimated sigma(s)" against PEtab's count of estimated parameters.
* **Inputs driven by condition parameters** (`u5 := piecewise(0, time - PdBu_time < 0, PdBu_dose)` with a flag, a dose and a time per
  condition): make the inputs plain parameters (0) that events set (`at (time >= T_PREEQUIL + 24): u5 = 1`) and put the conditions
  in `Events.py`. The six condition columns vanish from the model. Time restarts at the switch in PEtab, so the data times are
  shifted by the pre-equilibration length.
* **Replicate measurements** at one time are several rows of that (condition, time) in the data file (`emit.write_condition_data_csv`);
  a blank is "not measured" and the engine skips it.
* **A reaction that is always off** (`pu2 * u2` with `u2 := 0`) is dropped, with its parameter.
* **Pre-equilibration by integration** agrees with PEtab's steady-state solve at the nominal point and at most random points in the
  bounds, but this model is bistable and stiff in places: a few parameter sets make even the reference fail (CVODE
  "too much work") or land on another fixed point. Compare at several random sets and say which ones fail.

## 8. Pitfalls found on Sneyd_PNAS2002 (one shared sigma, a clean problem)

* **A power in a PEtab formula is `^`**, which Python reads as XOR: `reference._eval` translates it to `**`. Write it as `^` in an
  Antimony rule (`open_probability := (0.9 * A + 0.1 * O)^4`).
* **One noise parameter for every measurement across conditions** (`noiseParameters = sigma` on all rows): every loss element carries
  the same `"sigma_block"`, giving one profiled sigma for the whole problem (a per-replicate sigma, the default, is the wrong model).
* **A compartment with no species** (`default`) needs no volume parameter: pass `exclude=["V_Default"]` to `emit.write_parameters_csv`.
* **The benchmark's `simulations.tsv` can be wrong**: Sneyd's has two conditions swapped. Compare it by key (condition, time), not by
  row position (the tables are in different orders), and settle a disagreement with an independent implementation (AMICI's PEtab importer
  at the nominal parameters), not by editing the reference.
* **Time a warm AMICI run**: the first fit of a model compiles it (minutes); `benchmarks/` timings must come from a warm cache.

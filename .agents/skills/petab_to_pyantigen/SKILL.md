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

## 9. Pitfalls found on Smith_BMCSystBiol2013 (a 367-reaction network that is mostly one pattern)

* **Find the pattern before writing anything.** Normalise every species name in the reactions (a regex over the modification bits)
  and group the reactions by what they change; the `R<n>_<name>` rate parameters of one group are then a base value times a factor
  per modification (here 150 FOXO1 rate constants are twelve numbers: a base rate and 0.1, 0.5, 10, 22, 3 per modification). One
  loop plus a handful of named factors replaces them, and `tests/test_smith.py` checks every reaction (below).
* **A rate law in amounts needs the volume back.** PyAntiGen rewrites every species in a rate expression as a concentration,
  `X` -> `(X / V_Compartment)`, for every rate type (also `custom`). An SBML whose laws are amounts (substanceOnly species, a
  bimolecular `k * A * B`) therefore writes each species as `V_Compartment * X`; the module's `amounts()` does it in one place. The
  generated Antimony is then `V * (X / V)`: exact, and ugly. Compartment sizes enter this model once (R37r), as a ratio.
* **Species that no reaction changes are parameters.** Boundary species (`Ins`, `extracellular_ROS`, `Mt`), enzymes held constant
  (`PP2A`, `SCF`, `Proteasome`, `SGK`) and `E2F1` are parameters with the SBML's initial amount; a condition sets them in
  `Update_parameters.py` like any other parameter. Species and rules the observables and rates never use are dropped.
* **Events after the last datum are not simulated.** The SBML's insulin pulses end at `t_ins` (>= every measurement of the
  condition) and at 2880; `Events.py` is empty and `test_no_event_can_change_a_measured_value` says why.
* **The reference simulator had three gaps** (fixed in `reference.py`): `r.reset()` puts a constant global parameter back to its SBML
  value, so a condition's `k4` or `indicator_jnk` set before the reset was lost (`_reapply` sets them again after it);
  a species override must be an INITIAL value (`init(E2F1)`), and `resetToOrigin()` (not `resetAll()`) undoes it for the next
  condition; a species in an observable formula is its amount when the SBML says `hasOnlySubstanceUnits`. PEtab parameters that are
  not in the SBML (`sc_*`) are available to the observable formulas. Check the reference against `simulations.tsv` first (Smith's
  matches to 9e-8, so it is not swapped like Sneyd's).
* **A blank cell in conditions.tsv means "the model's own value"** (E2F1 is 300 where the table is blank, 150 elsewhere):
  `emit.effective_condition_values` resolves it, `emit.varying` finds the columns that differ, and a column constant over all conditions
  (`k4`, `kminus4`, `k_irs1_basal_syn`) belongs in the parameters file.
* **`noiseFormula = 1.0`** with no `noiseParameters` column and no noise parameter in the table: every sigma is declared
  (`sigma_method: "fixed"`), so nothing is profiled and LM works. `compare_fits.score` and `amici_direct` score that as a declared sigma.
* **Test the right-hand sides, not only the observables.** `test_every_rate_equals_the_sbml_at_random_states` loads the generated
  Antimony into RoadRunner and compares all 92 derivatives with the SBML's at random amounts (1e-11): a wrong factor in a rarely
  observed state shows up there and in no observable.
* **`python -m petab2pyantigen.verify --petab ...`**, not `python -m petab2pyantigen verify --petab ...` (argparse rejects the option
  after the sub-command).
* **AMICI and a `piecewise` threshold** (`piecewise(x - t, x > t, 0)`, the kinetics of Akt and PKC): AMICI integrates it as an event
  and CVODE sometimes fails where the state crosses the threshold ("root after reinitialization", status -99), whether it does is
  knife-edge in the parameters (a 1e-16 change decides it) and worse at tight tolerances. The fix is a looser solver tolerance tried in turn:
  `analytic_gradient: {"rtol": 1e-8, "atol": 1e-10, "retry": [[1e-6, 1e-8], [1e-5, 1e-7], [1e-4, 1e-6]]}`; `amici_direct.py` takes
  `--rtol --atol --retry`. Rewriting it as `(x + abs(x)) / 2` is exact but much worse: AMICI hangs on the sensitivities (the second
  derivative of |x| is a delta). A smooth `(x + sqrt(x^2 + eps^2)) / 2` (eps = 0.1 molecule against a threshold of 200: observables
  within 1e-4, objective within 0.4) integrates everywhere and made L-BFGS-B and LM converge from starts where the exact model's failed
  (one reached the lowest objective found); it is an option (`SMITH_PIP3_SWITCH_EPS`), not PEtab's model, so score its fits on the
  exact model (`fit_problem.py` labels the route "smoothed switch"). Look for "the analytic gradient failed N time(s)" in the fit log.
* **A gradient that is exact can still stall L-BFGS-B**: from one start PyAntiGen's L-BFGS-B ended after zero iterations (the line search
  steps to a corner of the box, objective 4e18, and then cannot recover) while its gradient equalled AMICI's exact one, and scipy's own
  L-BFGS-B on the same gradient converged. LM is the robust route here. Check a stalled fit against `amici_direct.py` before
  blaming the model.
* **LM can be discarded at the noise floor**: at an objective of ~2.4e4 the engine's default projected-gradient tolerance (0.01 nats)
  is below the noise of the integration, the steps are rejected, and the fit is reported "not converged" and not cached.
  `options: {"pg_tol": 2.0}` (1e-4 of the objective) fixes that; a fit that fails because every trial step failed to integrate is a real failure.
* **The reference can be the weak link**: RoadRunner cannot integrate the best LM parameters (CVODE too much work), AMICI can.
  `compare_fits` says so; check the number with AMICI's own likelihood (PEtab importer) instead of discarding the fit.
* **Budget the slow routes**: a RoadRunner evaluation of 36 conditions is ~1.4 s (Nelder-Mead budget 3,000 here, not 20,000), an AMICI
  gradient ~6 s; give each route a stated budget (`--max-evals`) and report "budget" where it ends there.

## 10. Pitfalls found on Schwen_PONE2014 (log10 observables, priors, initial amounts of fitted parameters)

* **`observableTransformation = log10`**: write the observable as a log10 RULE of the model (`observable_IR1 := log10(scale * (IR1 + IR1in + offset))`)
  and store the LOG10 of every measurement in the data file; the engine compares the two on that scale, which is PEtab's (normal noise on the
  transformed value). PEtab's density also carries the Jacobian `1 / (y ln 10)`, a constant of the data: `reference.log_likelihood` and
  `compare_fits.score` include it in -log L, the engine's NLL does not, and nothing that is fitted sees it. `verify.compare` compares an
  observable on its own scale (it transforms PEtab's simulation) and takes `column_of=` to pick the column from the whole measurement row.
  `tests/test_schwen.py::test_the_data_are_the_log10_of_petabs_measurements` checks that 10 ** column is PEtab's measurement.
* **One PEtab observable with per-series observable parameters** (`observable_Insulin` with `offset_nExpID1;scaleElisa_nExpID1;fragments;km_nExpID1`
  in one block of conditions and `..._nExpID2...` in another): one rule per series (`observable_Insulin_nExpID<k>`), each over its own
  parameters. Nothing is renamed at run time and every parameter is fitted by its PEtab name.
* **Two noise parameters, each shared by several observables AND several conditions** (`IR_obs_std` for three observables in seven
  conditions, `std` for the ELISA series in twelve): `sigma_pool` on the observable entries and `sigma_block` on the loss elements,
  together. Check the engine's line against PEtab's count: `2 block(s), 286 points, k=28 parameters + 2 estimated sigma(s) = 30`.
* **PEtab bounds an estimated noise parameter** (`IR_obs_std` <= 0.056); the engine's profiled sigma has no bound. At a poor start the ML
  sigma can be outside (0.109): `compare_fits.score` puts it on the bound (`sigma_at_bound`), `score(..., bounded=False)` is the engine's
  objective. At the optimum it is inside here, so the two agree.
* **Objective priors** (`parameterScaleNormal`, mean and sd in decades, on six rate constants): the engine has them only on the AMICI route,
  `analytic_gradient["priors"] = [{"param", "ref": 10**mean, "sigma": sd, "side": "both"}]`, so RoadRunner's Nelder-Mead optimizes the data
  term alone (they are small: sd = 3 decades). AMICI's own likelihood has none either: `amici_direct.py` adds `0.5 z^2` and its gradient.
  `reference.prior_terms` is PEtab's term. When checking a gradient against finite differences, hook `AmiciGradient.linear_gradient`,
  not `LossGradient.gradient`: the priors are added after the latter.
* **An initial amount that is an expression of fitted parameters** (`Rec1 = ini_R1`, `Rec2 = ini_R1 * ini_R2fold`): write the expression
  in `_InitialConditions.csv` (`Rec1_ExtracellularMedium,ini_R1`), so AMICI gets the initial assignment and differentiates it, and set the
  species again in `Update_opt_parameters` for RoadRunner (the engine resets the model BEFORE it sets the fitted values). A condition's own
  amount (`init_Ins`) is set on the species in `Update_parameters`; AMICI gets it as a fixed initial amount, which is right.
* **PyAntiGen bug found here (fixed in `Scaled_amici.py`)**: the AMICI runner read the initial state sensitivities from
  `model.get_initial_state_sensitivities()`, which in AMICI 1.1.0 puts the entries of every listed parameter in the FIRST parameter's slot
  once several are listed: the gradient with respect to `ini_R1` and `ini_R2fold` was exactly zero and the first parameter's was wrong.
  It now reads them from a run that stops at t0 (`tests/engine/test_amici_initial_sensitivity.py` in PyAntiGen fails without the fix).
  Any problem whose fitted parameters feed initial amounts was affected; check a new gradient against finite differences of the reference.
* **A condition measured only at t = 0** (the first dose): RoadRunner's `simulate` needs two times, and the engine a window longer than a
  point: the reference adds a second time, `solver_settings` simulates to one minute.
* **`parameter_scale` as a dict** (`{name: "lin" | "log10"}`) for a mix; `amici_direct.py` handles both (it used to refuse a `lin` parameter).
* **Species nothing observes and nothing reads** (`Uptake1`, `Uptake2`, running totals) are dropped with their reactions.
* **A small model is not made smaller by a loop**: here PyAntiGen's definition is 2.2x smaller than SBML + conditions in bytes but 3.1x LARGER
  than clean Tellurium Antimony; the saving is in the measurement table (36 kB of repeated ids against an 11 kB data file).
* **Where an AMICI gradient call goes (Schwen, 19 conditions, 0.43 s)**: AMICI 0.41 s, loss 5 ms, replay Jacobian 18 ms. Before `ParameterReplay.jacobian`
  had its selection-matrix check it was 0.25 s (57 RoadRunner replays per model per call) and a project that re-read its data CSV inside
  `Loss_config` per call cost 0.3 s more: cache anything a `loss_config` function reads (`functools.lru_cache`). Profile with `cProfile` around
  `AmiciGradient.linear_gradient` after one warm call; time a fit's s/eval only from a warm AMICI cache.

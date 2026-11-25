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
* **`six`** is needed by Tellurium but no longer declared by it (it is in requirements.txt).

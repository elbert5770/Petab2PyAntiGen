"""Write the data files of a PyAntiGen project from a PEtab problem.

The reactions are code (a module of f-string loops) and are written by hand or by
an agent; everything that is *data* in PEtab becomes one of three small files:

  antimony_models/<Model>/<Model>_parameters.csv         nominal values
  antimony_models/<Model>/<Model>_InitialConditions.csv   SBML initial amounts
  data/<Model>.csv                                        measurements, one row per time

PyAntiGen conventions applied here:
  * a species is ``<species>_<Compartment>`` (``K27me0K36me0_Default``); compartment
    names are CamelCase with no underscore;
  * every compartment has an explicit volume parameter ``V_<Compartment>``;
  * species are substance-only, so initial values are amounts (concentration x volume).
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

import libsbml
import pandas as pd

from .petab_reader import PetabProblem
from .reference import _reapply, _set, nominal_values


def camel_case(name: str) -> str:
    """``default`` -> ``Default``; ``brain_isf`` -> ``BrainIsf`` (no underscores allowed)."""
    return "".join(p[:1].upper() + p[1:] for p in re.split(r"[_\W]+", name) if p)


def species_name(species_id: str, compartment_id: str) -> str:
    return f"{species_id}_{camel_case(compartment_id)}"


class _Sbml:
    """Keeps the SBMLDocument alive: a libsbml Model dangles once its document is collected."""

    def __init__(self, problem: PetabProblem):
        self.doc = libsbml.readSBMLFromFile(str(problem.sbml_path))
        self.model = self.doc.getModel()


def _model(problem: PetabProblem) -> libsbml.Model:
    h = _Sbml(problem)
    h.model._keepalive = h          # the model carries its document with it
    return h.model


def assignment_targets(problem: PetabProblem) -> set[str]:
    """Ids set by SBML rules: they are rules in the Antimony model, not parameters."""
    return {r.getVariable() for r in _model(problem).getListOfRules()}


def preequilibration_values(problem: PetabProblem) -> dict[str, float]:
    """Condition overrides that hold while the model equilibrates (the model starts there)."""
    if "preequilibrationConditionId" not in problem.measurements:      # Smith: every simulation starts from the SBML state
        return {}
    pre = problem.measurements["preequilibrationConditionId"].dropna().unique()
    out: dict[str, float] = {}
    for cid in pre:
        for col, val in problem.conditions.loc[cid].items():
            if col not in ("conditionId", "conditionName") and pd.notna(val):
                out[col] = float(val)
    return out


def write_parameters_csv(problem: PetabProblem, path: str | Path, comments: dict[str, str] | None = None,
                         extra: dict[str, float] | None = None, exclude=()) -> list[str]:
    """Nominal values of every SBML parameter, plus one ``V_<Compartment>`` per compartment.

    PEtab's nominalValue wins over the value stored in the SBML (as in PEtab). ``exclude`` names SBML parameters the
    PyAntiGen model does not have (a condition's switches and doses, which become events; a rate that is always zero);
    ``extra`` adds parameters the SBML does not have (observable scale factors, inputs an event sets).
    Returns the parameter names written.
    """
    m = _model(problem)
    nominal = nominal_values(problem)
    nominal.update(preequilibration_values(problem))
    skip = assignment_targets(problem) | set(exclude)
    rows = []
    for c in m.getListOfCompartments():
        if f"V_{camel_case(c.getId())}" in set(exclude):          # a compartment without species (Sneyd's `default`)
            continue
        rows.append((f"V_{camel_case(c.getId())}", c.getSize(), "volume of compartment " + c.getId()))
    for p in m.getListOfParameters():
        if p.getId() in skip:
            continue
        rows.append((p.getId(), nominal.get(p.getId(), p.getValue()), ""))
    for k, v in (extra or {}).items():
        rows.append((k, v, ""))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Parameter", "Value", "Units", "Comment"])
        for name, val, comment in rows:
            w.writerow([name, repr(float(val)), "", (comments or {}).get(name, comment)])
    return [r[0] for r in rows]


def write_initial_conditions_csv(problem: PetabProblem, path: str | Path) -> list[str]:
    """SBML initial amounts, keyed by the PyAntiGen species name."""
    m = _model(problem)
    comp_size = {c.getId(): c.getSize() for c in m.getListOfCompartments()}
    rows = []
    for s in m.getListOfSpecies():
        amount = s.getInitialAmount() if s.isSetInitialAmount() else s.getInitialConcentration() * comp_size[s.getCompartment()]
        rows.append((species_name(s.getId(), s.getCompartment()), amount))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Species", "InitialCondition", "Units", "Comment"])
        for name, amount in rows:
            w.writerow([name, repr(float(amount)), "", ""])
    return [r[0] for r in rows]


def write_wide_data_csv(problem: PetabProblem, path: str | Path) -> pd.DataFrame:
    """Measurements pivoted to one row per time and one column per observable.

    PEtab's long table repeats ``observableId``, both condition ids and the noise
    parameter on every row; one condition and one noise parameter per problem
    leaves only (time x observable) numbers.  Requires a single (preequilibration,
    simulation) condition pair and no replicate (time, observable) duplicates.
    """
    meas = problem.measurements
    pairs = meas[["preequilibrationConditionId", "simulationConditionId"]].drop_duplicates()
    if len(pairs) != 1:
        raise NotImplementedError("write_wide_data_csv needs exactly one condition pair")
    wide = meas.pivot(index="time", columns="observableId", values="measurement")
    wide.insert(0, "time", wide.index)
    wide.to_csv(path, index=False, float_format="%.12g")
    return wide


# --------------------------------------------------------------------------------------------------------------------
# conditions
# --------------------------------------------------------------------------------------------------------------------

def effective_condition_values(problem: PetabProblem) -> pd.DataFrame:
    """What each PEtab condition does to the model: one row per condition, one column per quantity.

    For every condition the SBML is loaded, parameters.tsv nominal values are applied, the condition's overrides are set
    and the model is reset, so SBML rules and initial assignments that depend on the overrides are evaluated (Zhao's six
    0/1 "trigger" parameters select the parameters and initial infections of each region and stage this way). The row holds
    every global parameter and every species' initial amount.  Columns that differ between conditions are the ones a
    PyAntiGen replicate has to set (see :func:`varying`).
    """
    import roadrunner
    nominal = nominal_values(problem)
    r = roadrunner.RoadRunner(str(problem.sbml_path))
    params = list(r.model.getGlobalParameterIds())
    species = list(r.model.getFloatingSpeciesIds()) + list(r.model.getBoundarySpeciesIds())
    rows = {}
    cond_ids = list(dict.fromkeys(problem.measurements["simulationConditionId"]))
    for cid in cond_ids:
        r.resetToOrigin()                   # not resetAll(): it keeps a species override (init(E2F1)) of the previous condition
        for pid in set(params) & set(nominal):
            try:
                r[pid] = nominal[pid]
            except RuntimeError:                     # an assignment-rule target
                pass
        row = problem.conditions.loc[cid]
        cond = {col: float(val) for col, val in row.items() if col not in ("conditionId", "conditionName") and pd.notna(val)}
        for col, val in cond.items():
            _set(r, col, val)
        r.reset()
        _reapply(r, set(params), nominal, cond)
        rows[cid] = {**{p: float(r[p]) for p in params}, **{f"{s}(0)": float(r[s]) for s in species}}
    return pd.DataFrame.from_dict(rows, orient="index")


def varying(table: pd.DataFrame, rtol: float = 1e-12) -> list[str]:
    """Columns of :func:`effective_condition_values` whose value differs between conditions."""
    return [c for c in table.columns
            if not all(abs(v - table[c].iloc[0]) <= rtol * max(abs(table[c].iloc[0]), 1.0) for v in table[c])]


def write_condition_table_csv(problem: PetabProblem, path: str | Path, columns: dict[str, str],
                              labels: dict[str, dict[str, str]] | None = None) -> pd.DataFrame:
    """The compact per-condition table a PyAntiGen ``Experiment.py`` loops over.

    ``columns`` maps a column of :func:`effective_condition_values` to the name written (``{"Total_Pop": "Total_Pop",
    "Unquarantined_Infected(0)": "U0"}``); ``labels`` adds columns derived from the condition id, as
    ``{"region": {"model1_Wuhan_Stage_I": "Wuhan", ...}}``.
    """
    eff = effective_condition_values(problem)
    out = eff[list(columns)].rename(columns=columns)
    for name, mapping in (labels or {}).items():
        out.insert(0, name, [mapping[c] for c in out.index])
    out.insert(0, "condition", out.index)
    out.to_csv(path, index=False, float_format="%.12g")
    return out


def write_long_data_csv(problem: PetabProblem, path: str | Path, observable: str, column: str,
                        condition_column: str = "condition") -> pd.DataFrame:
    """One observable's measurements as rows of (condition, time, value): the data of several conditions in one file."""
    m = problem.measurements
    m = m[m["observableId"] == observable]
    out = pd.DataFrame({condition_column: m["simulationConditionId"].values, "time": m["time"].values,
                        column: m["measurement"].values})
    out.to_csv(path, index=False, float_format="%.12g")
    return out


def write_condition_data_csv(problem: PetabProblem, path: str | Path, condition_column: str = "condition") -> pd.DataFrame:
    """Every measurement of every simulation condition, one row per (condition, time, replicate), one column per observable.

    PEtab's long table repeats ``observableId``, both condition ids and the observable and noise parameters on every row.
    Here an observable's measurement sits in its own column; the second, third... measurement of the same observable at
    the same time goes in a further row of that (condition, time), so replicates are kept. A blank is "not measured".
    """
    m = problem.measurements.copy()
    m["_rep"] = m.groupby(["simulationConditionId", "observableId", "time"]).cumcount()
    wide = m.pivot_table(index=["simulationConditionId", "time", "_rep"], columns="observableId", values="measurement",
                         aggfunc="first").reset_index().drop(columns="_rep")
    wide = wide.rename(columns={"simulationConditionId": condition_column})
    wide.columns.name = None
    wide.to_csv(path, index=False, float_format="%.12g")
    return wide

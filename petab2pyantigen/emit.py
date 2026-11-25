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
from .reference import nominal_values


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
    pre = problem.measurements["preequilibrationConditionId"].dropna().unique()
    out: dict[str, float] = {}
    for cid in pre:
        for col, val in problem.conditions.loc[cid].items():
            if col not in ("conditionId", "conditionName") and pd.notna(val):
                out[col] = float(val)
    return out


def write_parameters_csv(problem: PetabProblem, path: str | Path, comments: dict[str, str] | None = None,
                         extra: dict[str, float] | None = None) -> list[str]:
    """Nominal values of every SBML parameter, plus one ``V_<Compartment>`` per compartment.

    PEtab's nominalValue wins over the value stored in the SBML (as in PEtab).
    Returns the parameter names written.
    """
    m = _model(problem)
    nominal = nominal_values(problem)
    nominal.update(preequilibration_values(problem))
    skip = assignment_targets(problem)
    rows = []
    for c in m.getListOfCompartments():
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

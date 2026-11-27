"""Simulate a PEtab v1 problem the way the PEtab specification says, with RoadRunner.

This is the ground truth a PyAntiGen conversion is checked against: same SBML,
same parameters.tsv, same conditions (including pre-equilibration to steady
state), same observables.  It supports the common case -- one SBML model, plain
numeric or parameter-valued condition overrides, ``lin``/``log``/``log10``
observables with Gaussian noise -- and refuses (NotImplementedError) anything it
cannot score faithfully, instead of silently giving a wrong number.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
import libsbml
import roadrunner

from .petab_reader import PetabProblem

_MATH_NS = {k: getattr(np, k) for k in
            ("exp", "log", "log10", "sqrt", "sin", "cos", "tan", "abs", "pow")}
_MATH_NS.update(ln=np.log, pi=math.pi)


def _is_number(x) -> bool:
    try:
        float(x)
        return not (isinstance(x, float) and math.isnan(x))
    except (TypeError, ValueError):
        return False


def nominal_values(problem: PetabProblem) -> dict[str, float]:
    """parameters.tsv nominalValue, on the linear scale (as PEtab stores it)."""
    return {pid: float(v) for pid, v in problem.parameters["nominalValue"].items()}


def _condition_values(problem: PetabProblem, cond_id: str, pars: dict[str, float]) -> dict[str, float]:
    """Resolve one conditions.tsv row to {model entity: number}."""
    row = problem.conditions.loc[cond_id]
    out = {}
    for col, val in row.items():
        if col in ("conditionId", "conditionName") or pd.isna(val):
            continue
        out[col] = float(val) if _is_number(val) else pars[str(val)]
    return out


def _set(r: roadrunner.RoadRunner, target: str, value: float) -> None:
    """Apply one condition override. A species is given its INITIAL amount: ``r[species] = x`` followed by ``r.reset()`` is
    undone by the reset (it restores a floating species' initial value), so the override would silently be lost."""
    if target in r.model.getFloatingSpeciesIds() or target in r.model.getBoundarySpeciesIds():
        r[f"init({target})"] = value
    r[target] = value


def _reapply(r: roadrunner.RoadRunner, model_params: set, pars: dict, condition: dict) -> None:
    """Set the global parameters again after ``r.reset()``.

    RoadRunner's reset puts a (constant) global parameter back to its SBML value (Smith's ``k4``, ``indicator_jnk``), so what was
    set before it is lost. It is set before the reset as well, because SBML initial assignments may read it."""
    for pid in model_params & set(pars):
        try:
            r[pid] = pars[pid]
        except RuntimeError:
            pass
    for tgt, val in condition.items():
        if tgt in model_params:
            try:
                r[tgt] = val
            except RuntimeError:
                pass


def _eval(formula: str, ns: dict):
    # PEtab formulas (muParser syntax) write a power as ``^``, which Python reads as XOR
    return eval(formula.replace("^", "**"), {"__builtins__": {}}, {**_MATH_NS, **ns})


def simulate_petab(problem: PetabProblem, parameters: dict[str, float] | None = None,
                   abs_tol: float = 1e-14, rel_tol: float = 1e-12) -> pd.DataFrame:
    """Return the PEtab ``simulations`` table (measurements with a ``simulation`` column).

    ``parameters`` overrides nominalValue for the named parameters (linear scale).
    """
    pars = nominal_values(problem)
    pars.update(parameters or {})
    meas = problem.measurements

    r = roadrunner.RoadRunner(str(problem.sbml_path))
    r.integrator.absolute_tolerance = abs_tol
    r.integrator.relative_tolerance = rel_tol
    r.integrator.maximum_num_steps = 200000
    model_params = set(r.model.getGlobalParameterIds())
    species = list(r.model.getFloatingSpeciesIds()) + list(r.model.getBoundarySpeciesIds())
    doc = libsbml.readSBMLFromFile(str(problem.sbml_path))
    substance_only = {sp.getId() for sp in doc.getModel().getListOfSpecies() if sp.getHasOnlySubstanceUnits()}
    for pid in model_params & set(pars):
        # assignment-rule targets (e.g. inflow) are not settable; their drivers are
        try:
            r[pid] = pars[pid]
        except RuntimeError:
            pass

    rows = []
    group_cols = ["preequilibrationConditionId", "simulationConditionId"]
    meas = meas.copy()
    meas["preequilibrationConditionId"] = meas.get("preequilibrationConditionId", pd.Series(dtype=object))
    for (pre, sim), block in meas.groupby(group_cols, dropna=False, sort=False):
        r.resetToOrigin()                   # not resetAll(): it keeps a species override (init(E2F1)) of the previous condition
        for pid in model_params & set(pars):
            try:
                r[pid] = pars[pid]
            except RuntimeError:
                pass
        # The model starts in its FIRST condition: the pre-equilibration condition if there is one, else the simulation
        # condition. Its values are applied before the initial state is evaluated, because SBML initial assignments may
        # depend on them (Zhao's initial infections depend on the stage and region "triggers").
        has_pre = isinstance(pre, str) and bool(pre)
        first = _condition_values(problem, pre if has_pre else sim, pars)
        for tgt, val in first.items():
            _set(r, tgt, val)
        r.reset()
        _reapply(r, model_params, pars, first)
        if has_pre:
            r.steadyState()
            for tgt, val in _condition_values(problem, sim, pars).items():
                _set(r, tgt, val)
        times = np.array(sorted(block["time"].unique()), dtype=float)
        if times[0] != 0.0:
            times = np.concatenate([[0.0], times])
        if len(times) < 2:                  # RoadRunner needs a start and one more time (a condition measured only at t = 0)
            times = np.concatenate([times, [times[0] + 1.0]])
        # a species in a formula is its concentration, unless the SBML says hasOnlySubstanceUnits (then its amount)
        res = r.simulate(times=times, selections=["time"] + [s if s in substance_only else f"[{s}]" for s in species])
        sim_df = pd.DataFrame(res, columns=["time"] + species)
        glob = {p: r[p] for p in r.model.getGlobalParameterIds()}
        for _, m in block.iterrows():
            obs = problem.observables.loc[m["observableId"]]
            line = sim_df.loc[np.isclose(sim_df["time"], m["time"])].iloc[0]
            ns = {**pars, **glob, **{s: line[s] for s in species}, **placeholders(m, pars, "observable", "observableParameters")}
            y = float(_eval(obs["observableFormula"], ns))
            rows.append({**m.to_dict(), "simulation": y})
    return pd.DataFrame(rows)


def placeholders(row, pars: dict, kind: str, column: str) -> dict:
    """``observableParameter<i>_<obs>`` / ``noiseParameter<i>_<obs>`` of one measurement row: each name in the row's
    semicolon-separated ``column`` is a number or the id of a parameter (looked up in ``pars``)."""
    raw = row.get(column)
    if raw is None or (isinstance(raw, float) and math.isnan(raw)) or str(raw) == "":
        return {}
    values = str(raw).split(";")
    return {f"{kind}Parameter{i}_{row['observableId']}": float(v) if _is_number(v) else pars[v.strip()]
            for i, v in enumerate(values, 1)}


def _transform(y, kind: str):
    if kind == "lin":
        return y
    return {"log": np.log, "log10": np.log10}[kind](y)


def log_likelihood(problem: PetabProblem, sims: pd.DataFrame,
                   parameters: dict[str, float] | None = None) -> float:
    """PEtab log-likelihood of ``sims`` (Gaussian noise, noiseFormula evaluated)."""
    pars = nominal_values(problem)
    pars.update(parameters or {})
    total = 0.0
    for _, row in sims.iterrows():
        obs = problem.observables.loc[row["observableId"]]
        if str(obs.get("noiseDistribution", "normal")) != "normal":
            raise NotImplementedError("only normal noise is supported")
        ns = {**pars, **placeholders(row, pars, "noise", "noiseParameters")}
        sigma = float(_eval(str(obs["noiseFormula"]), ns))
        kind = str(obs.get("observableTransformation", "lin"))
        res = _transform(row["measurement"], kind) - _transform(row["simulation"], kind)
        total += -0.5 * (math.log(2 * math.pi * sigma ** 2) + (res / sigma) ** 2)
        # PEtab's density is of the measurement, not of its logarithm: the Jacobian 1 / (y ln b) is a constant of the data
        if kind == "log":
            total -= math.log(row["measurement"])
        elif kind == "log10":
            total -= math.log(row["measurement"] * math.log(10))
    return total


def prior_terms(problem: PetabProblem, parameters: dict[str, float] | None = None) -> tuple[float, float]:
    """``(nll, constant)`` of PEtab's objective priors at ``parameters``: ``0.5 ((x - mean) / sd)^2`` for each estimated parameter whose
    ``objectivePriorType`` is ``parameterScaleNormal`` (x is the parameter on its own scale, log10 for most), and the normal density's
    constant ``ln(sd sqrt(2 pi))``. A parameter with no prior contributes nothing; another prior type is not implemented."""
    pars = nominal_values(problem)
    pars.update(parameters or {})
    table = problem.parameters
    nll = const = 0.0
    if "objectivePriorType" not in table:
        return nll, const
    for pid, row in table.iterrows():
        kind = row.get("objectivePriorType")
        if not isinstance(kind, str) or not kind or int(row["estimate"]) != 1:
            continue
        if kind != "parameterScaleNormal":
            raise NotImplementedError(f"objective prior {kind!r} on {pid}")
        mean, sd = (float(v) for v in str(row["objectivePriorParameters"]).split(";"))
        x = {"log10": math.log10, "log": math.log, "lin": lambda v: v}[str(row.get("parameterScale", "lin"))](pars[pid])
        nll += 0.5 * ((x - mean) / sd) ** 2
        const += math.log(sd * math.sqrt(2 * math.pi))
    return nll, const

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
import re

import numpy as np
import pandas as pd
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
    r[target] = value


def _eval(formula: str, ns: dict):
    return eval(formula, {"__builtins__": {}}, {**_MATH_NS, **ns})


def simulate_petab(problem: PetabProblem, parameters: dict[str, float] | None = None,
                   abs_tol: float = 1e-14, rel_tol: float = 1e-12) -> pd.DataFrame:
    """Return the PEtab ``simulations`` table (measurements with a ``simulation`` column).

    ``parameters`` overrides nominalValue for the named parameters (linear scale).
    """
    pars = nominal_values(problem)
    pars.update(parameters or {})
    meas = problem.measurements
    if "observableParameters" in meas and meas["observableParameters"].notna().any():
        raise NotImplementedError("observableParameters placeholders are not supported")

    r = roadrunner.RoadRunner(str(problem.sbml_path))
    r.integrator.absolute_tolerance = abs_tol
    r.integrator.relative_tolerance = rel_tol
    r.integrator.maximum_num_steps = 200000
    model_params = set(r.model.getGlobalParameterIds())
    species = list(r.model.getFloatingSpeciesIds()) + list(r.model.getBoundarySpeciesIds())
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
        r.resetAll()
        for pid in model_params & set(pars):
            try:
                r[pid] = pars[pid]
            except RuntimeError:
                pass
        # The model starts in its FIRST condition: the pre-equilibration condition if there is one, else the simulation
        # condition. Its values are applied before the initial state is evaluated, because SBML initial assignments may
        # depend on them (Zhao's initial infections depend on the stage and region "triggers").
        has_pre = isinstance(pre, str) and bool(pre)
        for tgt, val in _condition_values(problem, pre if has_pre else sim, pars).items():
            _set(r, tgt, val)
        r.reset()
        if has_pre:
            r.steadyState()
            for tgt, val in _condition_values(problem, sim, pars).items():
                _set(r, tgt, val)
        times = np.array(sorted(block["time"].unique()), dtype=float)
        if times[0] != 0.0:
            times = np.concatenate([[0.0], times])
        res = r.simulate(times=times, selections=["time"] + [f"[{s}]" for s in species])
        sim_df = pd.DataFrame(res, columns=["time"] + species)
        glob = {p: r[p] for p in r.model.getGlobalParameterIds()}
        for _, m in block.iterrows():
            obs = problem.observables.loc[m["observableId"]]
            line = sim_df.loc[np.isclose(sim_df["time"], m["time"])].iloc[0]
            ns = {**glob, **{s: line[s] for s in species}}
            y = float(_eval(obs["observableFormula"], ns))
            rows.append({**m.to_dict(), "simulation": y})
    return pd.DataFrame(rows)


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
        formula = str(obs["noiseFormula"])
        placeholders = re.findall(r"noiseParameter(\d+)_" + re.escape(row["observableId"]), formula)
        ns = dict(pars)
        if placeholders:
            vals = str(row["noiseParameters"]).split(";")
            for i, v in enumerate(vals, 1):
                ns[f"noiseParameter{i}_{row['observableId']}"] = float(v) if _is_number(v) else pars[v]
        sigma = float(_eval(formula, ns))
        kind = str(obs.get("observableTransformation", "lin"))
        res = _transform(row["measurement"], kind) - _transform(row["simulation"], kind)
        total += -0.5 * (math.log(2 * math.pi * sigma ** 2) + (res / sigma) ** 2)
    return total

"""Schwen_PONE2014: the PyAntiGen project reproduces the PEtab problem (insulin receptors in hepatocytes; observables on a log10 scale,
two estimated noise parameters, six log-normal priors, four ELISA series with their own observable parameters)."""
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))

from petab2pyantigen import emit
from petab2pyantigen.petab_reader import load_petab
from petab2pyantigen.reference import log_likelihood, nominal_values, prior_terms, simulate_petab
from petab2pyantigen.verify import compare

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT / "problems" / "Schwen_PONE2014"
MODEL = "Schwen_PONE2014"


def _series(row):
    """The PyAntiGen column of a PEtab measurement: Insulin is one observable in PEtab and four (one per ELISA series) here."""
    k = re.search(r"nExpID(\d)", str(row.get("observableParameters") or ""))
    return f"{row['observableId']}_nExpID{k.group(1)}" if k else row["observableId"]


def test_generated_model():
    """12 reactions (PEtab's 14, less the two Uptake running totals nobody observes), 9 species in one compartment, 7 rules; the two
    receptor types come from one loop."""
    subprocess.run([sys.executable, "Model_generate.py"], cwd=PROJECT, check=True, capture_output=True)
    text = (ROOT / "antimony_models" / MODEL / f"{MODEL}_reactions.txt").read_text()
    assert len(re.findall(r"^\w+ : .*;", text, flags=re.M)) == 12
    assert len(re.findall(r"^substanceOnly species \w+ in ExtracellularMedium$", text, flags=re.M)) == 9
    rules = (ROOT / "antimony_models" / MODEL / f"{MODEL}_rules.txt").read_text().splitlines()
    assert len(rules) == 7 and all(":= log10(" in r for r in rules)
    assert "Binding2 : Ins_ExtracellularMedium + Rec2_ExtracellularMedium -> IR2_ExtracellularMedium; ka1 * ka2fold" in text


def test_petab_reference_matches_the_benchmarks_own_simulations(schwen_petab):
    """All 286 rows agree (a condition measured only at t = 0 needs the reference to simulate past it)."""
    problem = load_petab(schwen_petab)
    sims = simulate_petab(problem)
    key = ["observableId", "simulationConditionId", "time", "observableParameters"]
    ref = problem.simulations.drop_duplicates(key)[key + ["simulation"]]
    both = sims.drop_duplicates(key).merge(ref, on=key, suffixes=("", "_ref"))
    assert len(sims) == 286 and len(both) == len(sims.drop_duplicates(key))
    assert ((both["simulation"] - both["simulation_ref"]).abs() / both["simulation_ref"].abs()).max() < 1e-6


def test_every_rate_equals_the_sbml_at_random_states(schwen_petab):
    """The nine time derivatives of the generated model equal the SBML's at random amounts and random parameters."""
    import libsbml
    import roadrunner
    import tellurium as te
    from pyantigen.generate.antimony_utils import load_antimony_files

    cwd = os.getcwd()
    os.chdir(ROOT)
    try:
        subprocess.run([sys.executable, "Model_generate.py"], cwd=PROJECT, check=True, capture_output=True)
        mine = te.loada(load_antimony_files(MODEL, str(ROOT)))
    finally:
        os.chdir(cwd)
    problem = load_petab(schwen_petab)
    ref = roadrunner.RoadRunner(str(problem.sbml_path))
    doc = libsbml.readSBMLFromFile(str(problem.sbml_path))                # keep the document: a Model dangles without it
    comp = {s.getId(): s.getCompartment() for s in doc.getModel().getListOfSpecies()}
    name = lambda sid: f"{sid}_{emit.camel_case(comp[sid])}"                      # noqa: E731
    rng = np.random.default_rng(0)
    species = list(mine.model.getFloatingSpeciesIds())
    assert len(species) == 9
    ref_ids = list(ref.model.getFloatingSpeciesIds())
    by_name = {name(s): s for s in ref_ids}
    shared = set(ref.model.getGlobalParameterIds()) & set(mine.model.getGlobalParameterIds())
    for _ in range(3):
        for p, v in nominal_values(problem).items():
            if p in shared:
                v = float(v * 10 ** rng.uniform(-0.5, 0.5))
                ref[p] = v
                mine[p] = v
        for s in species:
            v = float(10 ** rng.uniform(0, 3))
            mine[s] = v
            ref[by_name[s]] = v
        want = dict(zip(ref_ids, ref.getRatesOfChange()))
        got = dict(zip(species, mine.getRatesOfChange()))
        for s in species:
            assert got[s] == pytest.approx(want[by_name[s]], rel=1e-9, abs=1e-12), s


def test_pyantigen_matches_petab_at_nominal_and_at_random_parameters(schwen_petab):
    """Each log10 observable equals PEtab's (on the log10 scale), at the nominal parameters and at random ones. The receptor numbers
    are initial amounts that are set from the fitted ini_R1 and ini_R2fold: a random set moves them too."""
    r = compare(schwen_petab, PROJECT, "Schwen", column_of=_series)
    assert r["n"] == 286
    assert r["max_abs"] < 1e-7
    assert r["llh_pyantigen"] == pytest.approx(r["llh_petab"], rel=1e-8)
    problem = load_petab(schwen_petab)
    nominal, table = nominal_values(problem), problem.parameters
    rng = np.random.default_rng(1)
    pars = {k: float(np.clip(nominal[k] * 10 ** rng.uniform(-0.5, 0.5), float(table.loc[k, "lowerBound"]), float(table.loc[k, "upperBound"])))
            for k in table.index if table.loc[k, "estimate"] == 1 and k not in ("IR_obs_std", "std")}
    r = compare(schwen_petab, PROJECT, "Schwen", pars, column_of=_series)
    assert r["max_abs"] < 1e-7
    assert abs(r["llh_petab"] - r["llh_pyantigen"]) < 1e-6 * abs(r["llh_petab"])


def test_the_data_are_the_log10_of_petabs_measurements(schwen_petab):
    """286 values, one per PEtab row; 10 ** column is the measurement (the Insulin series are four columns)."""
    problem = load_petab(schwen_petab)
    data = pd.read_csv(ROOT / "data" / f"{MODEL}.csv")
    measured = data.drop(columns=["condition", "time"])
    assert int(measured.notna().sum().sum()) == len(problem.measurements) == 286
    mine = sorted(float(10 ** v) for v in measured.to_numpy().ravel() if not math.isnan(v))
    theirs = sorted(problem.measurements["measurement"].astype(float))
    assert mine == pytest.approx(theirs, rel=1e-10)
    for _, row in problem.measurements.iterrows():
        sub = data[(data["condition"] == row["simulationConditionId"]) & np.isclose(data["time"], row["time"])]
        assert any(np.isclose(10 ** v, row["measurement"], rtol=1e-10) for v in sub[_series(row)].dropna())


def test_objective_at_nominal_is_petabs_with_two_profiled_sigmas(schwen_petab):
    """Two blocks (IR_obs_std for the cytometry observables of seven conditions, std for the ELISA of twelve), 286 points, 28 + 2
    parameters as PEtab. RoadRunner's objective has no priors: it equals the reference's NLL less the priors' term."""
    from compare_fits import score
    problem = load_petab(schwen_petab)
    sc = score(problem, nominal_values(problem))
    proc = subprocess.run([sys.executable, "Model_run.py", "--optimize", "Schwen", "--no-fit"], cwd=PROJECT,
                          capture_output=True, text=True, encoding="utf-8", errors="replace")
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "2 block(s), 286 points, k=28 parameters + 2 estimated sigma(s) = 30" in proc.stdout
    nll = float(re.search(r"concentrated NLL at optimum: (-?[\d.]+)", proc.stdout).group(1))
    assert nll == pytest.approx(sc["nll"] - sc["prior_nll"], abs=2e-3)       # the engine prints 3 decimals
    # at PEtab's nominal point (itself a fit) the profiled sigma is PEtab's own, and the -log L is PEtab's
    assert sc["sigma_at_bound"] == []
    petab_nll = -log_likelihood(problem, simulate_petab(problem))
    assert petab_nll == pytest.approx(sc["neg_log_l"] - sc["prior_nll"] - prior_terms(problem)[1], abs=0.02)


def test_the_priors_are_petabs(schwen_petab):
    """Six parameterScaleNormal priors, mean in decades and sd 3, as the AMICI route's log-normal priors; the reference's term is 0.5 z^2."""
    problem = load_petab(schwen_petab)
    sys.path.insert(0, str(PROJECT))
    import AntiGen_paths  # noqa: F401
    from Modules.Optimizer_settings import PRIORS
    table = problem.parameters
    petab_priors = {p: tuple(float(v) for v in str(table.loc[p, "objectivePriorParameters"]).split(";"))
                    for p in table.index if isinstance(table.loc[p, "objectivePriorType"], str)}
    assert len(petab_priors) == 6 and {i["param"] for i in PRIORS} == set(petab_priors)
    for item in PRIORS:
        mean, sd = petab_priors[item["param"]]
        assert math.log10(item["ref"]) == pytest.approx(mean, rel=1e-12) and item["sigma"] == sd and item["side"] == "both"
    nominal = nominal_values(problem)
    by_hand = sum(0.5 * ((math.log10(nominal[p]) - m) / sd) ** 2 for p, (m, sd) in petab_priors.items())
    assert prior_terms(problem)[0] == pytest.approx(by_hand, rel=1e-12)


def test_the_condition_table_is_the_insulin_each_experiment_starts_with(schwen_petab):
    problem = load_petab(schwen_petab)
    eff = emit.effective_condition_values(problem)
    table = pd.read_csv(ROOT / "data" / f"{MODEL}_conditions.csv").set_index("condition")
    assert len(table) == 19
    assert table["init_Ins"].to_dict() == pytest.approx(eff["init_Ins"].to_dict(), rel=1e-12)

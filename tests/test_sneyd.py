"""Sneyd_PNAS2002: the PyAntiGen project reproduces the PEtab problem (nine constant-stimulus conditions, one shared sigma)."""
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from petab2pyantigen.petab_reader import load_petab
from petab2pyantigen.reference import nominal_values, simulate_petab
from petab2pyantigen.verify import compare

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT / "problems" / "Sneyd_PNAS2002"


def test_generated_model():
    """Ten transitions between six states in one compartment, the volume explicit; L1, L3, L5 and the observable are rules."""
    subprocess.run([sys.executable, "Model_generate.py"], cwd=PROJECT, check=True, capture_output=True)
    text = (ROOT / "antimony_models" / "Sneyd_PNAS2002" / "Sneyd_PNAS2002_reactions.txt").read_text()
    assert len(re.findall(r"^\w+ : .*;", text, flags=re.M)) == 10
    assert len(re.findall(r"^substanceOnly species \w+_Membrane in Membrane$", text, flags=re.M)) == 6
    assert "compartment Membrane := V_Membrane" in text
    rules = (ROOT / "antimony_models" / "Sneyd_PNAS2002" / "Sneyd_PNAS2002_rules.txt").read_text()
    assert "L1 := k_1 / k1 / (l_2 / l2)" in rules
    assert "open_probability := (0.9 * IPR_A_Membrane + 0.1 * IPR_O_Membrane)^4" in rules


def test_petab_reference_matches_the_benchmarks_own_simulations(sneyd_petab):
    """Seven of the nine conditions agree with simulations.tsv to 1e-6. Conditions 5 and 6 (3 and 10 uM Ca2+) agree with the
    OTHER condition's rows there: the benchmark's simulations.tsv has the two swapped. Not a defect of the reference: AMICI's own
    PEtab importer gives the same log-likelihood (319.7918) and the same two curves as the reference, so the conversion is held to
    the SBML and conditions.tsv, not to those 30 rows."""
    problem = load_petab(sneyd_petab)
    sims = simulate_petab(problem)
    ref = problem.simulations
    key = ["observableId", "simulationConditionId", "time"]
    both = sims.merge(ref[key + ["simulation"]], on=key, suffixes=("", "_ref"))
    assert len(both) == len(sims) == 135
    swapped = {"Ca_dose_response__5", "Ca_dose_response__6"}
    ok = both[~both["simulationConditionId"].isin(swapped)]
    assert (ok["simulation"] - ok["simulation_ref"]).abs().max() < 1e-6
    five = both[both["simulationConditionId"] == "Ca_dose_response__5"].sort_values("time")["simulation"].values
    six_ref = both[both["simulationConditionId"] == "Ca_dose_response__6"].sort_values("time")["simulation_ref"].values
    assert np.abs(five - six_ref).max() < 1e-6


def test_pyantigen_matches_petab_at_nominal(sneyd_petab):
    r = compare(sneyd_petab, PROJECT, "Sneyd")
    assert r["n"] == 135
    assert r["max_rel"] < 1e-6
    assert r["llh_pyantigen"] == pytest.approx(r["llh_petab"], abs=1e-4)


def test_pyantigen_matches_petab_at_a_random_parameter_set(sneyd_petab):
    problem = load_petab(sneyd_petab)
    nominal, table = nominal_values(problem), problem.parameters
    rng = np.random.default_rng(1)
    pars = {k: float(np.clip(nominal[k] * 10 ** rng.uniform(-0.5, 0.5), float(table.loc[k, "lowerBound"]),
                             float(table.loc[k, "upperBound"])))
            for k in table.index if k != "sigma"}
    r = compare(sneyd_petab, PROJECT, "Sneyd", pars)
    assert r["max_rel"] < 1e-6


def test_objective_at_nominal_is_petabs_with_one_shared_sigma(sneyd_petab):
    """PEtab estimates one sigma for all 135 points: one profiled block, (n/2) ln(SSE/n), and PEtab's own parameter count, 15."""
    problem = load_petab(sneyd_petab)
    sims = simulate_petab(problem)
    n = len(sims)
    expected = 0.5 * n * math.log(((sims.measurement - sims.simulation) ** 2).sum() / n)
    proc = subprocess.run([sys.executable, "Model_run.py", "--optimize", "Sneyd", "--no-fit"], cwd=PROJECT,
                          capture_output=True, text=True, encoding="utf-8", errors="replace")
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "1 block(s), 135 points, k=14 parameters + 1 estimated sigma(s) = 15" in proc.stdout
    nll = float(re.search(r"concentrated NLL at optimum: (-?[\d.]+)", proc.stdout).group(1))
    assert nll == pytest.approx(expected, abs=1e-2)


def test_the_condition_table_is_what_the_petab_conditions_do_to_the_model(sneyd_petab):
    """Only the constant IP3 and Ca2+ differ between the nine conditions."""
    import pandas as pd
    from petab2pyantigen import emit
    problem = load_petab(sneyd_petab)
    eff = emit.effective_condition_values(problem)
    assert sorted(emit.varying(eff)) == ["Ca", "IP3"]
    table = pd.read_csv(ROOT / "data" / "Sneyd_PNAS2002_conditions.csv").set_index("condition")
    for column in ("IP3", "Ca"):
        assert table[column].to_dict() == pytest.approx(eff[column].to_dict(), rel=1e-12)

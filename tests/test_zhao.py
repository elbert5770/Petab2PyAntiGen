"""Zhao_QuantBiol2020: the PyAntiGen project reproduces the PEtab problem."""
import math
import re
import subprocess
import sys
from pathlib import Path

import pytest

from petab2pyantigen.petab_reader import load_petab
from petab2pyantigen.reference import nominal_values, simulate_petab
from petab2pyantigen.verify import compare

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT / "problems" / "Zhao_QuantBiol2020"
OBSERVABLE = {"observable_confirmed_infected": "Confirmed_Infected_Population"}


def test_generated_model():
    """Three reactions over four species in one compartment, the volume explicit; no triggers anywhere."""
    subprocess.run([sys.executable, "Model_generate.py"], cwd=PROJECT, check=True, capture_output=True)
    text = (ROOT / "antimony_models" / "Zhao_QuantBiol2020" / "Zhao_QuantBiol2020_reactions.txt").read_text()
    assert len(re.findall(r"^\w+ : .*;", text, flags=re.M)) == 3
    assert len(re.findall(r"^substanceOnly species \w+_Population in Population$", text, flags=re.M)) == 4
    assert "compartment Population := V_Population" in text
    assert "alpha * (Susceptible_Population/V_Population) * (Unquarantined_Infected_Population/V_Population) / Total_Pop" in text
    assert "Trigger" not in text


def test_petab_reference_matches_the_benchmarks_own_simulations(zhao_petab):
    """Unlike Zheng's, this benchmark's simulations.tsv is a real simulation, so it checks our reference."""
    problem = load_petab(zhao_petab)
    sims = simulate_petab(problem)
    ref = problem.simulations
    key = ["observableId", "simulationConditionId", "time"]
    both = sims.merge(ref[key + ["simulation"]], on=key, suffixes=("", "_ref"))
    assert len(both) == len(sims) == 82
    assert ((both["simulation"] - both["simulation_ref"]).abs() / both["simulation_ref"]).max() < 1e-6


def test_pyantigen_matches_petab_at_nominal(zhao_petab):
    r = compare(zhao_petab, PROJECT, "Zhao", observable_name=OBSERVABLE, strip_prefix="model1_")
    assert r["n"] == 82
    assert r["max_rel"] < 1e-7
    assert r["llh_pyantigen"] == pytest.approx(r["llh_petab"], abs=1e-4)


def test_fitted_parameters_reach_the_right_replicate(zhao_petab):
    """R_Stage_I_Wuhan means R in the Wuhan, Stage I replicate only: through Update_opt_parameters, not by name."""
    problem = load_petab(zhao_petab)
    pars = {p: v * 1.4 for p, v in nominal_values(problem).items() if p.startswith(("R_", "gamma_"))}
    r = compare(zhao_petab, PROJECT, "Zhao", pars, observable_name=OBSERVABLE, strip_prefix="model1_")
    assert r["max_rel"] < 1e-7
    unchanged = compare(zhao_petab, PROJECT, "Zhao", None, observable_name=OBSERVABLE, strip_prefix="model1_")
    assert abs(r["table"]["pyantigen"] - unchanged["table"]["pyantigen"]).max() > 1.0      # the parameters did act


def test_objective_at_nominal_is_petabs_with_one_sigma_per_dataset(zhao_petab):
    """PEtab estimates sd per dataset; PyAntiGen profiles sigma per replicate: sum of (n/2) ln(SSE/n) over the 7 datasets."""
    problem = load_petab(zhao_petab)
    sims = simulate_petab(problem)
    expected = sum(0.5 * len(g) * math.log(((g.measurement - g.simulation) ** 2).sum() / len(g))
                   for _c, g in sims.groupby("simulationConditionId"))
    proc = subprocess.run([sys.executable, "Model_run.py", "--optimize", "Zhao", "--no-fit"], cwd=PROJECT,
                          capture_output=True, text=True, encoding="utf-8", errors="replace")
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "7 block(s), 82 points" in proc.stdout
    nll = float(re.search(r"concentrated NLL at optimum: (-?[\d.]+)", proc.stdout).group(1))
    assert nll == pytest.approx(expected, abs=1e-2)


def test_the_condition_table_is_what_the_petab_conditions_do_to_the_model(zhao_petab):
    """The six trigger parameters, the rules and the initial assignments reduce to seven rows of independent values."""
    import pandas as pd
    from petab2pyantigen import emit
    problem = load_petab(zhao_petab)
    eff = emit.effective_condition_values(problem)
    assert {"R", "gamma_1", "gamma_2", "Total_Pop", "Unquarantined_Infected(0)", "Quarantined_Infected(0)",
            "Confirmed_Infected(0)"} <= set(emit.varying(eff))
    table = pd.read_csv(ROOT / "data" / "Zhao_QuantBiol2020_conditions.csv").set_index("condition")
    for column, name in {"R": "R", "gamma_1": "gamma_1", "gamma_2": "gamma_2", "Total_Pop": "Total_Pop",
                         "U0": "Unquarantined_Infected(0)", "Q0": "Quarantined_Infected(0)",
                         "C0": "Confirmed_Infected(0)"}.items():
        assert table[column].to_dict() == pytest.approx(eff[name].to_dict(), rel=1e-12), column

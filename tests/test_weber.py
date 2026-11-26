"""Weber_BMC2015: the PyAntiGen project reproduces the PEtab problem (observable scale parameters, fixed noise, time-switched inputs)."""
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
PROJECT = ROOT / "problems" / "Weber_BMC2015"


def test_generated_model():
    """Twenty reactions over seven species in one compartment, the volume explicit; no switches or doses in the model."""
    subprocess.run([sys.executable, "Model_generate.py"], cwd=PROJECT, check=True, capture_output=True)
    text = (ROOT / "antimony_models" / "Weber_BMC2015" / "Weber_BMC2015_reactions.txt").read_text()
    assert len(re.findall(r"^\w+ : .*;", text, flags=re.M)) == 20
    assert len(re.findall(r"^substanceOnly species \w+_Cyt in Cyt$", text, flags=re.M)) == 7
    assert "compartment Cyt := V_Cyt" in text
    assert "flag" not in text and "dose" not in text
    rules = (ROOT / "antimony_models" / "Weber_BMC2015" / "Weber_BMC2015_rules.txt").read_text()
    assert "yPKDpN0 := PKDDAGa_Cyt * scale_yPKDpN0" in rules


def test_petab_reference_matches_the_benchmarks_own_simulations(weber_petab):
    """Scale parameters in the observables and a pre-equilibration followed by clock-driven inputs: the reference agrees with
    the benchmark's simulations.tsv."""
    problem = load_petab(weber_petab)
    sims = simulate_petab(problem)
    ref = problem.simulations
    assert len(sims) == len(ref) == 135
    rel = np.abs(sims["simulation"].values - ref["simulation"].values) / np.abs(ref["simulation"].values)
    assert rel.max() < 1e-6


def test_pyantigen_matches_petab_at_nominal(weber_petab):
    r = compare(weber_petab, PROJECT, "Weber", strip_prefix="")
    assert r["n"] == 135
    assert r["max_rel"] < 1e-7
    assert r["llh_pyantigen"] == pytest.approx(r["llh_petab"], abs=1e-3)


def test_pyantigen_matches_petab_at_a_random_parameter_set(weber_petab):
    """The scale factors, the rates and the inputs' strengths all act: every estimated parameter moved by up to a factor 3."""
    problem = load_petab(weber_petab)
    nominal, table = nominal_values(problem), problem.parameters
    rng = np.random.default_rng(1)
    pars = {k: float(np.clip(nominal[k] * 10 ** rng.uniform(-0.5, 0.5), float(table.loc[k, "lowerBound"]),
                             float(table.loc[k, "upperBound"])))
            for k, e in table["estimate"].items() if e == 1 and not k.startswith("std")}
    r = compare(weber_petab, PROJECT, "Weber", pars, strip_prefix="")
    assert r["max_rel"] < 1e-6
    assert r["llh_pyantigen"] == pytest.approx(r["llh_petab"], rel=1e-6)


def test_objective_at_nominal_is_petabs_with_five_estimated_and_three_declared_sigmas(weber_petab):
    """Estimated noise parameters are profiled per observable (yPKDpN0 across both conditions); the three fixed ones
    (estimate = 0) are declared: n ln(sigma) + SSE / (2 sigma^2)."""
    problem = load_petab(weber_petab)
    sims = simulate_petab(problem)
    fixed = {"std_yPKDt": 1910527.835, "std_yPI4K3Bt": 887676.131, "std_yCERTt": 244727220.419}
    expected = 0.0
    for noise, g in sims.groupby(sims["noiseParameters"].astype(str)):
        sse, n = float(((g.measurement - g.simulation) ** 2).sum()), len(g)
        expected += (n * math.log(fixed[noise]) + sse / (2 * fixed[noise] ** 2)) if noise in fixed else 0.5 * n * math.log(sse / n)
    proc = subprocess.run([sys.executable, "Model_run.py", "--optimize", "Weber", "--no-fit"], cwd=PROJECT,
                          capture_output=True, text=True, encoding="utf-8", errors="replace")
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "8 block(s), 135 points, k=31 parameters + 5 estimated sigma(s) = 36" in proc.stdout     # PEtab estimates 36
    nll = float(re.search(r"concentrated NLL at optimum: (-?[\d.]+)", proc.stdout).group(1))
    assert nll == pytest.approx(expected, abs=1e-2)

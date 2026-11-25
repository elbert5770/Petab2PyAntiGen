"""Zheng_PNAS2012: the PyAntiGen project reproduces the PEtab problem."""
import re
import subprocess
import sys
from pathlib import Path

import pytest

from petab2pyantigen.reference import log_likelihood, simulate_petab
from petab2pyantigen.petab_reader import load_petab
from petab2pyantigen.verify import compare

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT / "problems" / "Zheng_PNAS2012"


def test_generated_reaction_network():
    """60 reactions (22 reversible methylation steps = 44, 1 inflow, 15 dilutions), 15 species, one compartment."""
    subprocess.run([sys.executable, "Model_generate.py"], cwd=PROJECT, check=True, capture_output=True)
    text = (ROOT / "antimony_models" / "Zheng_PNAS2012" / "Zheng_PNAS2012_reactions.txt").read_text()
    assert len(re.findall(r"^\w+ : .*;", text, flags=re.M)) == 60
    assert len(re.findall(r"^\w+_fwd : ", text, flags=re.M)) == 22
    assert len(re.findall(r"^\w+_rev : ", text, flags=re.M)) == 22
    assert len(re.findall(r"^substanceOnly species K27me\dK36me\d_Default in Default$", text, flags=re.M)) == 15
    assert "compartment Default := V_Default" in text
    assert "K27me3K36me3" not in text


def test_petab_reference_is_self_consistent(zheng_petab):
    """At parameters.tsv nominal values sigma is the maximum-likelihood estimate: chi2 == N."""
    problem = load_petab(zheng_petab)
    sims = simulate_petab(problem)
    sigma = float(problem.parameters.loc["sigma", "nominalValue"])
    chi2 = (((sims["measurement"] - sims["simulation"]) / sigma) ** 2).sum()
    assert chi2 == pytest.approx(len(sims), rel=1e-3)
    assert log_likelihood(problem, sims) == pytest.approx(278.3335, abs=1e-3)


def test_pyantigen_matches_petab_at_nominal(zheng_petab):
    r = compare(zheng_petab, PROJECT, "Zheng")
    assert r["n"] == 60
    assert r["max_abs"] < 1e-8
    assert r["llh_pyantigen"] == pytest.approx(r["llh_petab"], abs=1e-6)


def test_pyantigen_matches_petab_at_other_parameters(zheng_petab):
    problem = load_petab(zheng_petab)
    pars = {p: float(v) * 1.7 for p, v in problem.parameters["nominalValue"].items() if p != "sigma"}
    r = compare(zheng_petab, PROJECT, "Zheng", pars)
    assert r["max_abs"] < 1e-8


def test_objective_at_nominal_is_petabs_with_one_shared_sigma(zheng_petab):
    """sigma_pool: 15 observables, one block of 60 points, one sigma; NLL + n/2 (1 + ln 2 pi) == -llh(PEtab)."""
    import math
    problem = load_petab(zheng_petab)
    proc = subprocess.run([sys.executable, "Model_run.py", "--optimize", "Zheng", "--no-fit"], cwd=PROJECT,
                          capture_output=True, text=True, encoding="utf-8", errors="replace")
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "1 block(s), 60 points" in proc.stdout
    nll = float(re.search(r"concentrated NLL at optimum: (-?[\d.]+)", proc.stdout).group(1))
    sigma = float(re.search(r"= (0\.\d+)\s+\(n=60, estimated\)", proc.stdout).group(1))
    assert sigma == pytest.approx(float(problem.parameters.loc["sigma", "nominalValue"]), rel=1e-3)
    assert nll + 30 * (1 + math.log(2 * math.pi)) == pytest.approx(-log_likelihood(problem, simulate_petab(problem)), abs=1e-2)

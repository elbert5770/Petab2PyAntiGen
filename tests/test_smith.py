"""Smith_BMCSystBiol2013: the PyAntiGen project reproduces the PEtab problem (insulin signalling, ROS and FOXO1; 36 conditions,
declared sigmas, a 48-state FOXO1 network written as one loop)."""
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from petab2pyantigen import emit
from petab2pyantigen.petab_reader import load_petab
from petab2pyantigen.reference import nominal_values, simulate_petab
from petab2pyantigen.verify import compare

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT / "problems" / "Smith_BMCSystBiol2013"
MODEL = "Smith_BMCSystBiol2013"


def test_generated_model():
    """313 reactions (PEtab's 367: the two kinases of Pa and the 16 DNA-bound states' transcription are written once), 92 species
    in four compartments, the nine observables as rules; the FOXO1 network is 16 states in each of three places."""
    subprocess.run([sys.executable, "Model_generate.py"], cwd=PROJECT, check=True, capture_output=True)
    text = (ROOT / "antimony_models" / MODEL / f"{MODEL}_reactions.txt").read_text()
    assert len(re.findall(r"^\w+ : .*;", text, flags=re.M)) == 313
    assert len(re.findall(r"^substanceOnly species \w+ in \w+$", text, flags=re.M)) == 92
    assert len(re.findall(r"^substanceOnly species Foxo1_\w+ in \w+$", text, flags=re.M)) == 48
    rules = (ROOT / "antimony_models" / MODEL / f"{MODEL}_rules.txt").read_text().splitlines()
    assert len(rules) == 9
    assert "PTP_activ__2D := sc_PTP * PTP1B_Cytoplasm / (PTP1B_Cytoplasm + PTP1B_ox_Cytoplasm)" in rules


def test_petab_reference_matches_the_benchmarks_own_simulations(smith_petab):
    """Unlike Sneyd's, this simulations.tsv is a real simulation: all 62 rows agree (the SBML has events and boundary species,
    amounts and tiny compartment sizes; the reference gives a species' amount, not its concentration, when the SBML says so)."""
    problem = load_petab(smith_petab)
    sims = simulate_petab(problem)
    key = ["observableId", "simulationConditionId", "time"]
    both = sims.merge(problem.simulations[key + ["simulation"]], on=key, suffixes=("", "_ref"))
    assert len(both) == len(sims) == 62
    assert ((both["simulation"] - both["simulation_ref"]).abs() / both["simulation_ref"].abs()).max() < 1e-6


def test_every_rate_equals_the_sbml_at_random_states(smith_petab):
    """The 92 time derivatives of the generated model equal the SBML's at random amounts (all reactions, in particular the
    FOXO1 rates rebuilt from a base rate and a factor per modification), with the JNK/IKK and FOXO1 blocks switched on."""
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
    problem = load_petab(smith_petab)
    ref = roadrunner.RoadRunner(str(problem.sbml_path))
    doc = libsbml.readSBMLFromFile(str(problem.sbml_path))
    comp = {s.getId(): s.getCompartment() for s in doc.getModel().getListOfSpecies()}

    def name(sid):                                   # SBML species id -> PyAntiGen species name
        c = comp[sid]
        return f"{sid[len(c) + 1:] if sid.startswith(c + '_') else sid}_{emit.camel_case(c)}"

    by_name = {name(s): s for s in comp}
    ref_params = set(ref.model.getGlobalParameterIds())
    for p, v in nominal_values(problem).items():
        if p in ref_params:
            ref[p] = v
    for p, v in (("k4", 0.000333333), ("kminus4", 0.003), ("k_irs1_basal_syn", 130.0)):     # what every condition sets
        ref[p] = v
    for p, v in dict(Ins=5000.0, extracellular_ROS=25000.0, indicator_jnk=1.0, indicator_foxo=1.0, E2F1=150.0).items():
        ref[p] = v
        mine[p] = v
    species = list(mine.model.getFloatingSpeciesIds())
    assert len(species) == 92
    ref_ids = list(ref.model.getFloatingSpeciesIds())
    rng = np.random.default_rng(0)
    for _ in range(2):
        for s in species:
            v = float(rng.uniform(0.5, 1.5) * 10 ** rng.uniform(0, 4))
            mine[s] = v
            ref[by_name[s]] = v
        want = dict(zip(ref_ids, ref.getRatesOfChange()))
        got = dict(zip(species, mine.getRatesOfChange()))
        scale = max(abs(v) for v in want.values())
        for s in species:
            assert got[s] == pytest.approx(want[by_name[s]], rel=1e-8, abs=1e-9 * scale), s


def test_pyantigen_matches_petab_at_nominal(smith_petab):
    r = compare(smith_petab, PROJECT, "Smith")
    assert r["n"] == 62
    assert r["max_rel"] < 1e-6
    assert r["llh_pyantigen"] == pytest.approx(r["llh_petab"], rel=1e-8)


def test_pyantigen_matches_petab_at_a_random_parameter_set(smith_petab):
    problem = load_petab(smith_petab)
    nominal, table = nominal_values(problem), problem.parameters
    rng = np.random.default_rng(1)
    pars = {k: float(np.clip(nominal[k] * 10 ** rng.uniform(-0.5, 0.5), float(table.loc[k, "lowerBound"]),
                             float(table.loc[k, "upperBound"])))
            for k in table.index if table.loc[k, "estimate"] == 1}
    r = compare(smith_petab, PROJECT, "Smith", pars)
    assert r["max_rel"] < 1e-6


def test_objective_at_nominal_is_petabs_with_every_sigma_declared(smith_petab):
    """noiseFormula is 1.0 everywhere and no noise parameter is estimated: nothing is profiled (51 blocks, all declared), the
    loss is half the sum of squares, and the count is PEtab's own, 25 parameters."""
    problem = load_petab(smith_petab)
    sims = simulate_petab(problem)
    expected = 0.5 * ((sims.measurement - sims.simulation) ** 2).sum()
    proc = subprocess.run([sys.executable, "Model_run.py", "--optimize", "Smith", "--no-fit"], cwd=PROJECT,
                          capture_output=True, text=True, encoding="utf-8", errors="replace")
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "51 block(s), 62 points, k=25 parameters + 0 estimated sigma(s) = 25" in proc.stdout
    nll = float(re.search(r"concentrated NLL at optimum: (-?[\d.]+)", proc.stdout).group(1))
    assert nll == pytest.approx(expected, rel=1e-6)


def test_the_condition_table_is_what_the_petab_conditions_do_to_the_model(smith_petab):
    """Five quantities differ between the 36 conditions; the other columns of conditions.tsv (k4, kminus4, k_irs1_basal_syn)
    hold for every condition and are the values in the parameters file. E2F1 is 150 in the FOXO1 experiments and the model's
    300 elsewhere (a blank in conditions.tsv): a species override that a RoadRunner reset undoes unless it is an initial value."""
    problem = load_petab(smith_petab)
    eff = emit.effective_condition_values(problem)
    varying = set(emit.varying(eff))
    assert {"Ins(0)", "extracellular_ROS(0)", "indicator_jnk", "indicator_foxo", "E2F1(0)"} <= varying
    assert not varying & {"k4", "kminus4", "k_irs1_basal_syn"}
    table = pd.read_csv(ROOT / "data" / f"{MODEL}_conditions.csv").set_index("condition")
    for column, source in (("Ins", "Ins(0)"), ("extracellular_ROS", "extracellular_ROS(0)"), ("indicator_jnk", "indicator_jnk"),
                           ("indicator_foxo", "indicator_foxo"), ("E2F1", "E2F1(0)")):
        assert table[column].to_dict() == pytest.approx(eff[source].to_dict(), rel=1e-12)
    assert set(table["E2F1"]) == {150.0, 300.0}
    parameters = pd.read_csv(ROOT / "antimony_models" / MODEL / f"{MODEL}_parameters.csv").set_index("Parameter")["Value"]
    for p in ("k4", "kminus4", "k_irs1_basal_syn"):
        assert parameters[p] == pytest.approx(float(eff[p].iloc[0]), rel=1e-12)


def test_no_event_can_change_a_measured_value(smith_petab):
    """The SBML has an insulin pulse that ends at t_ins and a second one at 2880 min. Every condition's last measurement is at or
    before its t_ins, and 2880 is later than any measurement: the project has no events (Modules/Events.py) and still equals
    PEtab, whose simulation includes them (the tests above)."""
    import libsbml
    problem = load_petab(smith_petab)
    doc = libsbml.readSBMLFromFile(str(problem.sbml_path))          # keep the document: a Model dangles without it
    events = doc.getModel().getListOfEvents()
    assert [e.getId() for e in events] == ["insulin_stimulation_end", "insulin_restimulation_start", "insulin_restimulation_end"]
    last = problem.measurements.groupby("simulationConditionId")["time"].max()
    t_ins = problem.conditions["t_ins"].astype(float)
    assert (last <= t_ins.loc[last.index]).all()
    assert problem.measurements["time"].max() < 2880


def _generate(env):
    subprocess.run([sys.executable, "Model_generate.py"], cwd=PROJECT, check=True, capture_output=True,
                   env={**os.environ, "SMITH_PIP3_SWITCH_EPS": env})
    return (ROOT / "antimony_models" / MODEL / f"{MODEL}_reactions.txt").read_text()


def test_the_smoothed_switch_is_an_option_and_off_by_default(smith_petab):
    """SMITH_PIP3_SWITCH_EPS=0.1 writes the Akt/PKC switch as (x + sqrt(x^2 + eps^2)) / 2 instead of PEtab's piecewise: no AMICI event,
    observables within 1e-3 of PEtab's (it differs from max(x, 0) by at most eps / 2 = 0.05 molecules), and not equal to them."""
    try:
        exact = _generate("")
        assert "piecewise" in exact and "sqrt" not in exact
        smooth = _generate("0.1")                                # the files the simulation below uses
        assert "piecewise" not in smooth and smooth.count("sqrt(") == 2
        r = compare(smith_petab, PROJECT, "Smith")
        assert 1e-9 < r["max_rel"] < 1e-3
    finally:
        _generate("")                                            # leave the exact model in antimony_models/

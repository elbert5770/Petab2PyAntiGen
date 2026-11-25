"""How much smaller is the PyAntiGen project than the PEtab problem?

Files are grouped by what they are, because the groups compress very differently:

  PEtab                         PyAntiGen
  -----                         ---------
  model.xml                     antimony_modules/...        reactions and rules, as f-string loops
  observables.tsv               Modules/Loss_config.py      observables, noise
  conditions.tsv                Modules/Events.py + Solver_settings.py
  parameters.tsv                Modules/Optimizer_settings.py (bounds, scale) + <Model>_parameters.csv
  measurements.tsv              data/<Model>.csv            wide, one row per time
  problem.yaml, visualization   Modules/Experiment.py, Plots.py

A PyAntiGen project also carries fixed harness files (AntiGen_paths.py, Model_run.py ...)
that are the same for every model; they are listed separately so a one-model comparison
is not blamed for them.
"""
from __future__ import annotations

from pathlib import Path

from .petab_reader import PetabProblem

HARNESS = {"AntiGen_paths.py", "Model_run.py", "Experiment.py", "Observed_species.py",
           "Update_parameters.py", "Update_opt_parameters.py", "Model_generate.py", "__init__.py"}


def _size(paths) -> tuple[int, int]:
    """(bytes, non-blank lines)."""
    b = n = 0
    for p in paths:
        data = Path(p).read_bytes()
        b += len(data)
        n += sum(1 for ln in data.decode("utf-8", "replace").splitlines() if ln.strip())
    return b, n


def project_files(root: Path, model: str, modules: list[str]) -> dict[str, list[Path]]:
    """PyAntiGen source files of one problem, by group.  ``modules`` are paths under antimony_modules/."""
    prob = root / "problems" / model
    py = [p for p in prob.rglob("*.py") if "__pycache__" not in p.parts]
    am = root / "antimony_models" / model
    return {
        "reactions + rules (antimony_modules)": [root / "antimony_modules" / m for m in modules],
        "problem scripts (Modules/)": [p for p in py if p.name not in HARNESS],
        "parameters + initial values": [am / f"{model}_parameters.csv", am / f"{model}_InitialConditions.csv"],
        "measurements (data/)": [root / "data" / f"{model}.csv"],
        "harness (same for every model)": [p for p in py if p.name in HARNESS and p.name != "__init__.py"],
    }


def report(problem: PetabProblem, root: str | Path, model: str, modules: list[str],
           antimony_text: str | None = None) -> str:
    root = Path(root)
    lines = []
    petab = problem.files
    groups = {
        "SBML model": petab["sbml"],
        "parameters.tsv": petab["parameters"],
        "observables.tsv": petab["observables"],
        "conditions.tsv": petab["conditions"],
        "measurements.tsv": petab["measurements"],
        "problem.yaml + visualization": petab["yaml"] + petab["visualization"],
    }
    pb = pn = 0
    lines.append(f"PEtab problem {problem.name}")
    for name, paths in groups.items():
        b, n = _size(paths)
        pb, pn = pb + b, pn + n
        lines.append(f"  {name:34s} {b:9,d} bytes {n:6,d} lines")
    lines.append(f"  {'TOTAL':34s} {pb:9,d} bytes {pn:6,d} lines")
    if antimony_text is not None:
        lines.append(f"  (SBML exported to clean Antimony:  {len(antimony_text.encode()):9,d} bytes "
                     f"{sum(1 for l in antimony_text.splitlines() if l.strip()):6,d} lines)")

    lines.append(f"\nPyAntiGen project {model}")
    ab = an = hb = hn = 0
    for name, paths in project_files(root, model, modules).items():
        b, n = _size(paths)
        lines.append(f"  {name:34s} {b:9,d} bytes {n:6,d} lines")
        if name.startswith("harness"):
            hb, hn = b, n
        else:
            ab, an = ab + b, an + n
    lines.append(f"  {'TOTAL excluding harness':34s} {ab:9,d} bytes {an:6,d} lines   "
                 f"({100 * ab / pb:.0f}% of PEtab bytes, {100 * an / pn:.0f}% of lines)")
    lines.append(f"  {'TOTAL including harness':34s} {ab + hb:9,d} bytes {an + hn:6,d} lines")
    return "\n".join(lines)

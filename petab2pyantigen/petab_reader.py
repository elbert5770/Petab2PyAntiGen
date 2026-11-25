"""Load a PEtab v1 problem directory into plain pandas tables."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
import yaml


@dataclass
class PetabProblem:
    """The pieces of a PEtab v1 problem, as read from ``problem.yaml``."""

    directory: Path
    sbml_path: Path
    parameters: pd.DataFrame
    conditions: pd.DataFrame
    measurements: pd.DataFrame
    observables: pd.DataFrame
    simulations: pd.DataFrame | None = None
    files: dict = field(default_factory=dict)  # role -> list[Path]

    @property
    def name(self) -> str:
        return self.directory.parent.name if self.directory.name.startswith("v") else self.directory.name

    def file_sizes(self) -> dict[str, int]:
        """Bytes of every file that defines the problem (not the derived ones)."""
        return {p.name: p.stat().st_size for paths in self.files.values() for p in paths}


def _read(path: Path, index: str | None = None) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    return df.set_index(index, drop=False) if index else df


def load_petab(problem_yaml: str | Path) -> PetabProblem:
    """Read ``problem.yaml`` (PEtab format_version 1) and every table it names.

    ``problem_yaml`` may also be the directory that contains ``problem.yaml``.
    """
    problem_yaml = Path(problem_yaml)
    if problem_yaml.is_dir():
        problem_yaml = problem_yaml / "problem.yaml"
    root = problem_yaml.parent
    spec = yaml.safe_load(problem_yaml.read_text(encoding="utf-8"))
    if str(spec.get("format_version", "1")).split(".")[0] != "1":
        raise NotImplementedError("only PEtab format_version 1 problems are supported")
    if len(spec["problems"]) != 1:
        raise NotImplementedError("PEtab problems with more than one sub-problem")
    sub = spec["problems"][0]

    def paths(names):
        return [root / n for n in names]

    files = {
        "yaml": [problem_yaml],
        "parameters": paths([spec["parameter_file"]]),
        "conditions": paths(sub.get("condition_files", [])),
        "measurements": paths(sub.get("measurement_files", [])),
        "observables": paths(sub.get("observable_files", [])),
        "sbml": paths(sub.get("sbml_files", [])),
        "visualization": paths(sub.get("visualization_files", [])),
    }
    if len(files["sbml"]) != 1:
        raise NotImplementedError("PEtab problems with more than one SBML model")

    def concat(role, index):
        frames = [_read(p) for p in files[role]]
        df = pd.concat(frames, ignore_index=True)
        return df.set_index(index, drop=False) if index else df

    sims = root / "simulations.tsv"
    return PetabProblem(
        directory=root,
        sbml_path=files["sbml"][0],
        parameters=_read(files["parameters"][0], "parameterId"),
        conditions=concat("conditions", "conditionId"),
        measurements=concat("measurements", None),
        observables=concat("observables", "observableId"),
        simulations=_read(sims) if sims.exists() else None,
        files=files,
    )

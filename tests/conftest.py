import os
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
# PEtab benchmark collection; override with PETAB_BENCHMARK_DIR
BENCHMARK = Path(os.environ.get("PETAB_BENCHMARK_DIR", ROOT.parent / "petab-benchmark-problems"))


@pytest.fixture(scope="session")
def zheng_petab():
    d = BENCHMARK / "problems" / "Zheng_PNAS2012" / "v1"
    if not (d / "problem.yaml").exists():
        pytest.skip(f"PEtab Zheng_PNAS2012 not found at {d} (set PETAB_BENCHMARK_DIR)")
    return d


@pytest.fixture(scope="session")
def zhao_petab():
    d = BENCHMARK / "problems" / "Zhao_QuantBiol2020" / "v1"
    if not (d / "problem.yaml").exists():
        pytest.skip(f"PEtab Zhao_QuantBiol2020 not found at {d} (set PETAB_BENCHMARK_DIR)")
    return d

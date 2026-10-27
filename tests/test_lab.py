"""The browser lab engine (demo/lab_engine.js) must reproduce the papers and match Python."""
import json
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(shutil.which("node") is None, reason="node.js not installed")


@pytest.fixture(scope="module")
def lab():
    out = subprocess.run(["node", str(ROOT / "tests/lab_parity.js")], capture_output=True, text=True,
                         check=True).stdout
    return json.loads(out)


@pytest.mark.parametrize("case,time_min,energy_kwh,n", [
    ("none", 9.74, 0.1702, 1), ("tail", 5.50, 0.0539, 1), ("d45", 6.32, 0.0705, 1),
    ("head", 16.29, 0.4775, 3), ("spatial_tail", 5.95, 0.0635, 2),
    ("phase_none", 9.74, 0.1875, 1), ("phase_head", 16.29, 0.4950, 3)])
def test_trr_table_values_in_browser(lab, case, time_min, energy_kwh, n):
    r = lab["trr"][case]
    assert abs(r["T"] - time_min) < 0.01 and abs(r["E"] - energy_kwh) < 0.0006 and r["n"] == n


def test_trr_temporal_in_browser(lab):
    first, last = lab["trr"]["temporal"]
    assert 8.6 < first < 9.0 and 5.8 < last < 6.2


def test_sire_identical_to_python(lab):
    ref = json.loads((ROOT / "tests/data/sire_reference.json").read_text())
    assert lab["trb"]["same_decisions"]
    assert lab["trb"]["accepted"] == ref["kpis"]["accepted"]
    assert abs(lab["trb"]["gini"] - ref["kpis"]["gini"]) < 1e-9

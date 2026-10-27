"""Real-city layer: the browser engine (demo/geo_engine.js) must agree with Python."""
import json
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.skipif(shutil.which("node") is None, reason="node.js not installed")
@pytest.mark.skipif(not (ROOT / "tests/data/geo_reference.json").exists(),
                    reason="run python make_geo_demo.py first")
def test_js_engine_matches_python():
    out = subprocess.run(["node", str(ROOT / "tests/geo_parity.js")], capture_output=True,
                         text=True, check=True).stdout
    r = json.loads(out)
    assert r["max_node_wind_diff_mph"] < 1e-6
    assert r["sire_same"]
    for m, x in r["routes"].items():
        # identical route, or an exactly tied alternative under the planning wind
        assert x["same_nodes"] or x["dPlannedE_rel"] < 1e-9, m


def test_geo_corridor_avoids_buildings():
    import numpy as np
    from aamsim.geo import DALLAS, build_geo_corridor, load_buildings
    city = load_buildings(ROOT / "data/osm/dallas.osm.gz")
    gc = build_geo_corridor(city, DALLAS["start"], DALLAS["goal"])
    g, c = gc.grid, gc.corridor
    i = ((c.pos[:, 0] - g["x0"]) / g["res"]).astype(int)
    j = ((c.pos[:, 1] - g["y0"]) / g["res"]).astype(int)
    air = c.pos[:, 2] > 0
    assert np.all(c.pos[air, 2] > g["H"][i[air], j[air]])
    assert c.origins and c.destinations

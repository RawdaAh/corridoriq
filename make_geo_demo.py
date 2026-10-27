"""Build the real-city viewer (demo/city.html): OSM buildings of downtown Dallas,
the modular corridor laid over them, building-aware wind, wind-aware mission
planning with four planning winds, and SIRE corridor traffic - all computed live
in the browser by demo/geo_engine.js (the twin of the Python engine).

    python make_geo_demo.py
Also writes tests/data/geo_reference.json used by tests/test_geo.py to check
that the browser engine and the Python engine agree.
"""
from __future__ import annotations

import base64
import json
import math
from pathlib import Path

import numpy as np

from aamsim.corridor import Corridor
from aamsim.energy import geo_vehicle
from aamsim.geo import DALLAS, build_geo_corridor, load_buildings, load_roads, nearest_portal
from aamsim.geo_plan import run_geo_mission
from aamsim.sire import SireConfig, SireEngine, Trip
from aamsim.urbanwind import WindParams, solve_urban_wind

OSM = "data/osm/dallas.osm.gz"


def rounded_corridor(c: Corridor) -> Corridor:
    pos = np.round(c.pos, 2)
    return Corridor(c.spec, pos, c.ijk, c.arc_from, c.arc_to,
                    np.linalg.norm(pos[c.arc_to] - pos[c.arc_from], axis=1), c.arc_edge,
                    c.edge_nodes, c.origins, c.destinations)


def main():
    city = load_buildings(OSM, name=DALLAS["name"])
    gc = build_geo_corridor(city, DALLAS["start"], DALLAS["goal"])
    gc.corridor = rounded_corridor(gc.corridor)
    corr, g, fr = gc.corridor, gc.grid, gc.frame
    L, half = g["L"], g["half"]
    H16 = np.round(g["H"] * 10).astype("<u2")
    blds = []
    for b in city.buildings:
        xy = fr.to_corridor(b.xy)
        c = xy.mean(0)
        if -700 < c[0] < L + 700 and abs(c[1]) < 800:
            blds.append(dict(p=np.round(xy, 1).ravel().tolist(), h=round(b.height, 1),
                             b=round(b.base, 1), m=int(b.measured)))
    roads = []
    for en, major in load_roads(OSM, city.origin):
        xy = fr.to_corridor(en)
        if np.any((xy[:, 0] > -700) & (xy[:, 0] < L + 700) & (np.abs(xy[:, 1]) < 800)):
            roads.append(dict(p=np.round(xy, 1).ravel().tolist(), M=int(major)))
    o, d = nearest_portal(gc, "origin"), nearest_portal(gc, "destination")
    data = dict(
        meta=dict(city=DALLAS["name"], start=DALLAS["start_name"], goal=DALLAS["goal_name"],
                  start_ll=DALLAS["start"], goal_ll=DALLAS["goal"], length_m=round(L, 1),
                  bearing_deg=round(math.degrees(fr.bearing), 3), n_buildings=len(city.buildings),
                  measured_heights=int(sum(b.measured for b in city.buildings)),
                  over_ceiling=int(sum(b.height > 120 for b in city.buildings)),
                  removed_links=gc.removed_edges),
        grid=dict(x0=g["x0"], y0=g["y0"], nx=g["nx"], ny=g["ny"], res=g["res"],
                  H=base64.b64encode(H16.tobytes()).decode()),
        corridor=dict(pos=corr.pos.tolist(), edges=corr.edge_nodes.tolist(), ijk=corr.ijk.tolist(),
                      origins=corr.origins, destinations=corr.destinations, half=half,
                      z_levels=list(map(float, corr.spec.altitudes())), o=o, d=d,
                      module=[corr.spec.module_length, corr.spec.module_width]),
        vehicle=dict(cruise_mph=35.0, energy_j_per_ft=20.42, length_unit_ft=3.28084),
        wind=dict(ref_height=120.0, alpha=0.25, iterations=300, omega=1.85, wake_len=90.0,
                  wake_depth=0.65, wake_reach=24, max_speedup=1.8),
        buildings=blds, roads=roads)
    blob = json.dumps(data, separators=(",", ":"))
    engine = Path("demo/geo_engine.js").read_text()
    html = (Path("demo/city_template.html").read_text()
            .replace("/*__ENGINE__*/", engine).replace("/*__DATA__*/null", blob))
    Path("demo/city.html").write_text(html)
    print(f"demo/city.html: {len(html) / 1e6:.2f} MB, {len(blds)} buildings, {len(roads)} roads")

    # ---------------- parity reference for the JS engine
    gc.grid["H"] = H16.astype(float) / 10.0
    p = WindParams(speed_ms=8.0, direction_deg=30.0)
    uw = solve_urban_wind(gc.grid, corr.spec.altitudes(), p)
    rows, routes, _ = run_geo_mission(gc, p, o, d, uw=uw)
    from aamsim.urbanwind import wind_at_nodes
    wn = wind_at_nodes(uw, corr)
    from aamsim.energy import wind_link_costs
    veh = geo_vehicle()
    t_true = wind_link_costs(corr, veh, wn, "crab")[0][0]
    rng = np.random.default_rng(7)
    trips = [Trip(k + 1, int(rng.integers(1, 4)), int(rng.choice(corr.origins)),
                  int(rng.choice(corr.destinations)), int(rng.integers(0, 40))) for k in range(150)]
    eng = SireEngine(corr, veh, SireConfig(slot_sec=10.0, route_horizon=120, ground_transit=False),
                     link_time_sec=t_true)
    res = eng.run(trips)
    acc = res.trips.sort_values("trip_id").accepted.astype(int).tolist()
    ref = dict(params=dict(speed_ms=8.0, direction_deg=30.0), o=o, d=d,
               node_u=wn.u[0].tolist(), node_v=wn.v[0].tolist(),
               routes={m: dict(nodes=r["nodes"], E=r["e"][-1], T=r["t"][-1], PE=r["planned_e"][-1])
                       for m, r in routes.items()},
               trips=[t.__dict__ for t in trips], sire_accepted=acc)
    Path("tests/data").mkdir(parents=True, exist_ok=True)
    Path("tests/data/geo_reference.json").write_text(json.dumps(ref))
    Path("tests/data/geo_data.json").write_text(blob)
    print("tests/data/geo_reference.json written;", sum(acc), "of", len(acc), "SIRE trips accepted")


if __name__ == "__main__":
    main()

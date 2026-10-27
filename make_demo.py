"""Build the interactive replay demo (demo/index.html).

Runs small scenarios of all three studies plus the wind-information bridge,
exports geometry, routes, reservations and wind to JSON, and injects it into
``demo/template.html``. Open ``demo/index.html`` in a browser (it needs internet
access for the three.js script).

    python make_demo.py            # ~2-4 minutes
    python make_demo.py --fast     # fewer scenarios
"""
from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np

from aamsim.belief import METHODS, belief_corridor, run_mission, sample_missions
from aamsim.corridor import trb_sire_corridor, trr_wind_corridor
from aamsim.energy import VehicleModel, wind_link_costs
from aamsim.pareto import pareto_routes
from aamsim.rolling import PRIORITY, RHConfig, run_case
from aamsim.sire import SireConfig, generate_trips, run_sire
from aamsim.wind import trr_altitude_wind, trr_spatial_wind, trr_temporal_wind, zero_wind

R = lambda x, n=2: float(round(float(x), n))  # noqa: E731


def geom(corr):
    return dict(pos=[[R(a, 1) for a in p] for p in corr.pos],
                edges=corr.edge_nodes.tolist(), origins=corr.origins,
                destinations=corr.destinations, summary=corr.summary())


def wind_nodes(w, b=0):
    return dict(u=[R(x, 1) for x in w.u[b]], v=[R(x, 1) for x in w.v[b]])


def route_json(r):
    return dict(nodes=[int(n) for n in r.nodes], t=[R(x, 1) for x in r.times],
                e=[R(x / 3600, 2) for x in r.energies],  # Wh
                T=R(r.travel_min, 3), E=R(r.energy_kwh, 4))


def build_p1():
    corr = trr_wind_corridor()
    veh = VehicleModel()
    o, d = corr.node(0, 0, 0), corr.node(5, 0, 0)
    out = dict(corridor=geom(corr), origin=o, destination=d, cases={})
    for exp, maker in (("orientation", trr_altitude_wind), ("spatial", trr_spatial_wind)):
        rows = []
        for lab, w in [("No wind", zero_wind(corr))] + [(f"{a}°", maker(corr, a)) for a in range(0, 360, 45)]:
            rows.append(dict(label=lab, wind=wind_nodes(w),
                             routes=[route_json(r) for r in pareto_routes(corr, veh, w, o, d)]))
        out["cases"][exp] = rows
    temporal = []
    for dirn in (0, 180):
        for trend in ("calm_to_strong", "strong_to_calm"):
            w = trr_temporal_wind(corr, dirn, trend)
            costs = wind_link_costs(corr, veh, w)
            deps = []
            for k in range(10):
                R_ = pareto_routes(corr, veh, w, o, d, depart_sec=180 * k, costs=costs)
                deps.append(dict(dep_min=3 * k, routes=[route_json(r) for r in R_]))
            temporal.append(dict(label=f"{'Tailwind' if dirn == 0 else 'Headwind'}, {trend.replace('_', ' ')}",
                                 bin_sec=w.bin_sec,
                                 u=[[R(x, 1) for x in w.u[b]] for b in range(w.n_bins)],
                                 v=[[R(x, 1) for x in w.v[b]] for b in range(w.n_bins)],
                                 departures=deps))
    out["temporal"] = temporal
    return out


def build_p2(fast):
    corr = trb_sire_corridor()
    scen = [("Light demand", "400 trips, one direction", dict(demand=400)),
            ("Flow breakdown", "1,000 trips, one direction", dict(demand=1000)),
            ("Stress test", "2,000 trips, one direction", dict(demand=2000)),
            ("Two-way, shared links", "2,000 trips, 50/50 split, undirected key",
             dict(demand=2000, forward_share=0.5, rule="undirected")),
            ("Two-way, separated", "2,000 trips, 50/50 split, direction-separated key",
             dict(demand=2000, forward_share=0.5, rule="directed"))]
    if fast:
        scen = scen[:2]
    out = dict(corridor=geom(corr), slot_sec=60, scenarios=[])
    for lab, desc, p in scen:
        trips = generate_trips(corr, p["demand"], forward_share=p.get("forward_share", 1.0), seed=1001)
        res = run_sire(corr, trips, SireConfig(conflict_rule=p.get("rule", "undirected")))
        rows = []
        for r, legs in zip(res.trips.itertuples(), res.legs):
            nodes = None
            times = None
            if legs:
                nodes = [int(corr.arc_from[legs[0][0]])] + [int(corr.arc_to[a]) for a, _, _ in legs]
                times = [int(legs[0][1])] + [int(e) for _, _, e in legs]
            rows.append([int(r.operator), int(r.origin), int(r.destination), int(r.dispatch),
                         int(bool(r.accepted)), nodes, times,
                         R(r.travel_sec, 1) if r.accepted else None])
        k = res.kpis()
        out["scenarios"].append(dict(label=lab, desc=desc, kpis={a: (R(b, 4) if isinstance(b, float) else b)
                                                                  for a, b in k.items()}, trips=rows))
        print(f"  SIRE {lab}: SR={k['success_rate']:.3f}")
    return out


def build_p3(fast):
    cases = [("Base: 50 requests, peak", RHConfig())]
    if not fast:
        cases.append(("Wind: 15 mph headwind on forward flights", RHConfig(wind_speed_mph=15.0,
                                                                          wind_direction_deg=180.0)))
    out = dict(cases=[])
    for lab, cfg in cases:
        table, results, req, net = run_case(cfg)
        methods = {}
        for m, res in results.items():
            plans = []
            for rid, p in res["plans"].items():
                r = req.set_index("request_id").loc[rid]
                plans.append(dict(rid=int(rid), prio=r.priority, op=int(r.operator_id),
                                  path=[int(n) for n in p["node_path"]],
                                  times=[int(p["legs"][0][2])] + [int(x[3]) for x in p["legs"]],
                                  commit=int(p.get("commit_slot", 0)), delay=R(p["delay_min"], 2)))
            row = table[table.method == m].iloc[0]
            epochs = [dict(t=int(e["epoch"]), committed=int(e["committed"]),
                           tentative={str(k): [int(v[0])] for k, v in e.get("tentative", {}).items()})
                      for e in res["logs"]] if m != "offline" else []
            methods[m] = dict(plans=plans, epochs=epochs,
                              kpis=dict(served=int(row.served), service=R(row.service_rate, 3),
                                        weighted=R(row.weighted_service_rate, 3),
                                        mean_delay=R(row.mean_delay, 3), max_delay=R(row.max_delay, 2),
                                        emergency=R(row.emergency_service_rate, 3),
                                        replans=int(row.replans), max_wall=R(row.max_wall_sec, 2)))
        out["cases"].append(dict(
            label=lab, corridor=geom(net.corr), config=asdict(cfg),
            requests=[dict(rid=int(r.request_id), op=int(r.operator_id), o=int(r.origin),
                           d=int(r.destination), reveal=int(r.reveal), pref=int(r.preferred),
                           latest=int(r.latest), prio=r.priority) for r in req.itertuples()],
            methods=methods))
        print(f"  RH {lab}: done")
    out["priority_weights"] = PRIORITY
    return out


def build_p4(fast):
    corr = belief_corridor()
    veh = VehicleModel()
    ms = sample_missions(corr, 60, seed=0)
    picks = []
    for cls in ("headwind", "tailwind", "crosswind"):
        picks.append(next(m for m in ms if m.wind_class == cls))
    if fast:
        picks = picks[:1]
    out = dict(corridor=geom(corr), methods=METHODS, missions=[])
    for m in picks:
        rows, routes, truth, _ = run_mission(corr, veh, m, 15.0)
        out["missions"].append(dict(
            label=f"{m.wind_class.capitalize()} mission", wind_class=m.wind_class,
            wind_dir=R(m.wind_dir, 1), heading=m.heading, depart_sec=m.depart_bin * truth.bin_sec,
            bin_sec=truth.bin_sec,
            u=[[R(x, 1) for x in truth.u[b]] for b in range(truth.n_bins)],
            v=[[R(x, 1) for x in truth.v[b]] for b in range(truth.n_bins)],
            routes={k: dict(v, planned_t=[R(x, 1) for x in v["planned_t"]],
                            realized_t=[R(x, 1) for x in v["realized_t"]] if v["realized_t"] else None)
                    for k, v in routes.items()},
            results=[{k: (R(v, 4) if isinstance(v, float) else v) for k, v in r.items()} for r in rows]))
        print(f"  Belief {m.wind_class}: done")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true")
    ap.add_argument("--template", default="demo/template.html")
    ap.add_argument("--out", default="demo/index.html")
    ap.add_argument("--from-json", action="store_true", help="reuse demo/data.json, only rebuild HTML")
    a = ap.parse_args()
    t0 = time.time()
    if a.from_json:
        blob = Path("demo/data.json").read_text()
        Path(a.out).write_text(Path(a.template).read_text().replace("/*__DATA__*/null", blob))
        print(f"Rebuilt {a.out}")
        return
    data = dict(generated=time.strftime("%Y-%m-%d"), p1=build_p1(), p2=build_p2(a.fast),
                p3=build_p3(a.fast), p4=build_p4(a.fast))
    blob = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    Path("demo/data.json").write_text(blob)
    html = Path(a.template).read_text().replace("/*__DATA__*/null", blob)
    Path(a.out).write_text(html)
    print(f"Wrote {a.out} ({len(html) / 1e6:.1f} MB) in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()

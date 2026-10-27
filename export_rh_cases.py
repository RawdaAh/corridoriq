"""Run a set of rolling-horizon cases (Paper 3) with the Python ILP engine and export
complete schedules and epoch logs for the browser lab (demo/lab.html replays them).

    python export_rh_cases.py            # writes demo/rh_cases.json
"""
from __future__ import annotations

import json
import time
from dataclasses import asdict, replace
from pathlib import Path

from aamsim.rolling import PRIORITY, RHConfig, run_case

CASES = [
    ("base", "Base: 50 requests, central peak", {}),
    ("load_100", "Load: 100 requests", dict(n_requests=100)),
    ("freeze_4", "Freeze window 4 min", dict(freeze_min=4.0)),
    ("horizon_10", "Look-ahead 10 min", dict(horizon_min=10.0)),
    ("paths_1", "Single path (K = 1)", dict(k_paths=1)),
    ("two_peak", "Two demand peaks", dict(pattern="two_peak")),
    ("wind_head", "15 mph headwind on forward flights", dict(wind_speed_mph=15.0, wind_direction_deg=180.0)),
]


def export(key, label, over):
    cfg = replace(RHConfig(), **over)
    t0 = time.time()
    table, results, req, net = run_case(cfg)
    methods = {}
    lookup = req.set_index("request_id")
    for m, res in results.items():
        plans = []
        for rid, p in res["plans"].items():
            r = lookup.loc[rid]
            plans.append(dict(rid=int(rid), prio=r.priority, op=int(r.operator_id),
                              path=[int(n) for n in p["node_path"]],
                              times=[int(p["legs"][0][2])] + [int(x[3]) for x in p["legs"]],
                              commit=int(p.get("commit_slot", 0)), delay=round(float(p["delay_min"]), 3)))
        row = table[table.method == m].iloc[0]
        epochs = [dict(t=int(e["epoch"]), active=int(e["active"]), committed=int(e["committed"]),
                       new=int(e["new_commits"]), wall=round(float(e["wall_sec"]), 3),
                       tentative={str(k): int(v[0]) for k, v in e.get("tentative", {}).items()})
                  for e in res["logs"]] if m != "offline" else []
        methods[m] = dict(plans=plans, epochs=epochs, kpis={
            k: (None if row[k] != row[k] else round(float(row[k]), 4)) for k in
            ["served", "service_rate", "weighted_service_rate", "mean_delay", "p95_delay", "max_delay",
             "emergency_service_rate", "high_service_rate", "normal_service_rate", "replans",
             "max_wall_sec", "operator_rate_range"]})
    c = net.corr
    print(f"  {key}: {time.time() - t0:.0f} s", flush=True)
    return dict(key=key, label=label, config=asdict(cfg),
                corridor=dict(pos=c.pos.tolist(), edges=c.edge_nodes.tolist(),
                              origins=c.origins, destinations=c.destinations),
                requests=[dict(rid=int(r.request_id), op=int(r.operator_id), o=int(r.origin),
                               d=int(r.destination), reveal=int(r.reveal), pref=int(r.preferred),
                               latest=int(r.latest), prio=r.priority) for r in req.itertuples()],
                methods=methods)


def main():
    out = dict(priority_weights=PRIORITY, cases=[export(*c) for c in CASES])
    Path("demo/rh_cases.json").write_text(json.dumps(out, separators=(",", ":")))
    print("demo/rh_cases.json written")


if __name__ == "__main__":
    main()

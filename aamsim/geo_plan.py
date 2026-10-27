"""Energy-optimal planning on the real-city corridor and evaluation in the true wind.

Static wind snapshot (one bin): the energy-optimal route under a belief is a
shortest path on link energy (ties broken by time). This is the min-energy member
of the TRR Pareto set. The route is then flown in the ground-truth field.
"""
from __future__ import annotations

import heapq
import math

import numpy as np

from .energy import J_PER_KWH, VehicleModel, geo_vehicle, wind_link_costs
from .urbanwind import WindParams, geo_beliefs, solve_urban_wind

METHODS = ["ground_truth", "forecast", "profile", "wind_agnostic"]
LABELS = {"ground_truth": "Ground truth", "forecast": "Coarse forecast",
          "profile": "Profile wind", "wind_agnostic": "Wind-agnostic"}


def min_energy_route(corr, e_cost, t_cost, origin, destination):
    """Dijkstra on (energy, time) lexicographic. Returns arc list or None."""
    out = corr.out_arcs()
    best = {origin: (0.0, 0.0)}
    pred = {}
    pq = [(0.0, 0.0, origin)]
    ground = corr.ground
    while pq:
        e, t, u = heapq.heappop(pq)
        if (e, t) > best.get(u, (math.inf, math.inf)):
            continue
        if u == destination:
            break
        for a in out[u]:
            v = int(corr.arc_to[a])
            if ground[v] and v != destination:
                continue  # no intermediate landings
            ce, ct = e_cost[a], t_cost[a]
            if not np.isfinite(ce):
                continue
            cand = (e + ce, t + ct)
            if cand < best.get(v, (math.inf, math.inf)):
                best[v] = cand
                pred[v] = a
                heapq.heappush(pq, (cand[0], cand[1], v))
    if destination not in best:
        return None
    arcs, n = [], destination
    while n != origin:
        a = pred[n]
        arcs.append(int(a))
        n = int(corr.arc_from[a])
    return arcs[::-1]


def fly(corr, arcs, t_cost, e_cost, veh):
    t = np.concatenate([[0.0], np.cumsum(t_cost[arcs])])
    e = np.concatenate([[0.0], np.cumsum(e_cost[arcs])])
    return t, e


def classify(direction_deg, heading_deg):
    rel = abs((direction_deg - heading_deg + 180) % 360 - 180)
    return "tailwind" if rel <= 45 else ("headwind" if rel >= 135 else "crosswind")


def run_geo_mission(gc, params: WindParams, origin, destination, veh=None, model="paper",
                    uw=None):
    """Plan with the four beliefs, evaluate in the truth. Returns (rows, routes, uw)."""
    veh = veh or geo_vehicle()
    corr = gc.corridor
    if uw is None:
        uw = solve_urban_wind(gc.grid, corr.spec.altitudes(), params)
    beliefs = geo_beliefs(uw, corr, params)
    t_true, e_true = (c[0] for c in wind_link_costs(corr, veh, beliefs["ground_truth"], model))
    heading = 0.0 if corr.pos[destination, 0] > corr.pos[origin, 0] else 180.0
    rows, routes = [], {}
    for m in METHODS:
        tb, eb = (c[0] for c in wind_link_costs(corr, veh, beliefs[m], model))
        arcs = min_energy_route(corr, eb, tb, origin, destination)
        if arcs is None:
            rows.append(dict(method=m, feasible=False))
            continue
        tp, ep = fly(corr, arcs, tb, eb, veh)
        tr, er = fly(corr, arcs, t_true, e_true, veh)
        dist = float(corr.arc_len[arcs].sum())
        nodes = [int(corr.arc_from[arcs[0]])] + [int(corr.arc_to[a]) for a in arcs]
        z = corr.pos[nodes, 2]
        rows.append(dict(method=m, feasible=bool(np.isfinite(tr[-1])), distance_m=dist,
                         planned_s=tp[-1], realized_s=tr[-1], planned_wh=ep[-1] / 3600,
                         realized_wh=er[-1] / 3600, wh_per_km=er[-1] / 3600 / (dist / 1000),
                         mean_alt_m=float(np.average(0.5 * (z[1:] + z[:-1]), weights=np.diff(tr))),
                         ground_speed_ms=dist / tr[-1], arrival_error_s=tr[-1] - tp[-1]))
        routes[m] = dict(nodes=nodes, t=tr.tolist(), e=er.tolist(), planned_t=tp.tolist(),
                         planned_e=ep.tolist())
    ref = next((r for r in rows if r["method"] == "ground_truth" and r.get("feasible")), None)
    agn = next((r for r in rows if r["method"] == "wind_agnostic" and r.get("feasible")), None)
    for r in rows:
        r.update(direction_deg=params.direction_deg, speed_ms=params.speed_ms,
                 wind_class=classify(params.direction_deg, heading))
        if ref and r.get("feasible"):
            r["extra_energy_pct"] = 100 * (r["realized_wh"] / ref["realized_wh"] - 1)
        if agn and r.get("feasible"):
            r["saving_vs_agnostic_pct"] = 100 * (1 - r["realized_wh"] / agn["realized_wh"])
            if ref:
                full = agn["realized_wh"] - ref["realized_wh"]
                r["saving_recovered_pct"] = (100 * (agn["realized_wh"] - r["realized_wh"]) / full
                                             if full > 1e-9 else np.nan)
    return rows, routes, uw

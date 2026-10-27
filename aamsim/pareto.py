"""Time-dependent bi-objective (travel time, energy) label-correcting search.

TRR paper (Paper 1), Algorithm 1-2. Each label stores cumulative travel time,
cumulative energy and its predecessor. Labels are kept per (node, wind time
bin) and only non-dominated labels survive. Link costs depend on the wind bin
in which the link is entered. Labels exceeding the battery or the planning
horizon are discarded. Paths are elementary (no node revisits).

Equal-cost labels (same time and energy within a tolerance) are treated as
duplicates, so geometrically different but cost-identical routes are counted
once, as in the paper's tables.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass

import numpy as np

from .corridor import Corridor
from .energy import J_PER_KWH, VehicleModel, wind_link_costs

TOL_T = 1e-6
TOL_E = 1e-6


@dataclass
class ParetoRoute:
    nodes: list            # node sequence
    arcs: list             # arc ids
    depart_sec: float
    times: list            # cumulative time (s) at each node, relative to departure
    energies: list         # cumulative energy (J) at each node, incl. take-off
    travel_sec: float
    energy_j: float

    @property
    def travel_min(self):
        return self.travel_sec / 60.0

    @property
    def energy_kwh(self):
        return self.energy_j / J_PER_KWH


def _dominates(t1, e1, t2, e2):
    """(t1, e1) weakly dominates or duplicates (t2, e2)."""
    return t1 <= t2 + TOL_T and e1 <= e2 + TOL_E


def pareto_routes(corr: Corridor, veh: VehicleModel, wind, origin: int, destination: int,
                  depart_sec: float = 0.0, horizon_sec: float = 3600.0,
                  model: str = "paper", link_mask=None, costs=None, max_labels: int = 2_000_000,
                  ground_transit: bool = False):
    """All Pareto-optimal routes from origin to destination for one departure time.

    ``link_mask`` (bool, per arc) closes links (used for the congestion experiment).
    ``costs`` may pass precomputed (time, energy) arrays of shape (B, A).
    ``ground_transit=False`` forbids landing at an intermediate ground node and
    relaunching (only the origin and destination are ground nodes on a route).
    """
    if costs is None:
        costs = wind_link_costs(corr, veh, wind, model)
    tt, ee = costs
    out = corr.out_arcs()
    af, at = corr.arc_from, corr.arc_to
    open_arc = np.ones(corr.n_arcs, bool) if link_mask is None else np.asarray(link_mask, bool)
    e_to = veh.takeoff_energy_j()
    ground = corr.ground
    cap = veh.battery_kwh * J_PER_KWH

    # label arrays (grow as python lists for simplicity)
    L_t, L_e, L_node, L_pred, L_arc, L_vis, L_alive = [], [], [], [], [], [], []

    def new_label(t, e, node, pred, arc, vis):
        L_t.append(t); L_e.append(e); L_node.append(node); L_pred.append(pred)
        L_arc.append(arc); L_vis.append(vis); L_alive.append(True)
        return len(L_t) - 1

    buckets: dict = {}          # (node, bin) -> list of label ids
    dest_front: list = []       # non-dominated labels at the destination
    root = new_label(0.0, e_to, origin, -1, -1, 1 << origin)
    queue = deque([root])
    while queue:
        lid = queue.popleft()
        if not L_alive[lid]:
            continue
        u, t_u, e_u, vis = L_node[lid], L_t[lid], L_e[lid], L_vis[lid]
        b = int(wind.bin_at(depart_sec + t_u))
        for a in out[u]:
            if not open_arc[a]:
                continue
            v = int(at[a])
            if (vis >> v) & 1:
                continue
            if not ground_transit and ground[v] and v != destination:
                continue
            dt, de = tt[b, a], ee[b, a]
            if not np.isfinite(dt):
                continue
            t_v, e_v = t_u + dt, e_u + de
            if t_v > horizon_sec or e_v > cap:
                continue
            # bound against destination front
            if any(_dominates(L_t[d], L_e[d], t_v, e_v) for d in dest_front):
                continue
            if v == destination:
                dest_front = [d for d in dest_front
                              if not (t_v <= L_t[d] + TOL_T and e_v <= L_e[d] + TOL_E)]
                nid = new_label(t_v, e_v, v, lid, a, vis | (1 << v))
                dest_front.append(nid)
                continue
            key = (v, int(wind.bin_at(depart_sec + t_v)))
            bucket = buckets.setdefault(key, [])
            if any(L_alive[x] and _dominates(L_t[x], L_e[x], t_v, e_v) for x in bucket):
                continue
            for x in bucket:
                if L_alive[x] and t_v <= L_t[x] + TOL_T and e_v <= L_e[x] + TOL_E:
                    L_alive[x] = False
            bucket[:] = [x for x in bucket if L_alive[x]]
            nid = new_label(t_v, e_v, v, lid, a, vis | (1 << v))
            bucket.append(nid)
            queue.append(nid)
            if len(L_t) > max_labels:
                raise RuntimeError("label limit exceeded; refine the network or horizon")

    routes = []
    for d in sorted(dest_front, key=lambda i: (L_t[i], L_e[i])):
        arcs, nodes, times, energies = [], [], [], []
        x = d
        while x != -1:
            nodes.append(L_node[x]); times.append(L_t[x]); energies.append(L_e[x])
            if L_arc[x] >= 0:
                arcs.append(L_arc[x])
            x = L_pred[x]
        nodes.reverse(); times.reverse(); energies.reverse(); arcs.reverse()
        routes.append(ParetoRoute(nodes, arcs, depart_sec, times, energies, L_t[d], L_e[d]))
    return routes


def evaluate_route(corr: Corridor, veh: VehicleModel, wind, arcs, depart_sec: float,
                   model: str = "paper", costs=None):
    """Fly a fixed arc sequence through a (possibly different) wind field.

    Returns cumulative times and energies at every node; inf if a link is impassable.
    """
    if costs is None:
        costs = wind_link_costs(corr, veh, wind, model)
    tt, ee = costs
    t, e = 0.0, veh.takeoff_energy_j()
    times, energies = [t], [e]
    for a in arcs:
        b = int(wind.bin_at(depart_sec + t))
        t += tt[b, a]
        e += ee[b, a]
        times.append(t); energies.append(e)
    return np.array(times), np.array(energies)


# ----------------------------------------------------------------- route descriptors

def altitude_label(z0: float, z1: float) -> str:
    if z0 == z1:
        return f"{z0:g}"
    lo, hi = sorted((z0, z1))
    return f"{lo:g}-{hi:g}"


def route_profile(corr: Corridor, route: ParetoRoute, n_blocks: int = 5) -> dict:
    """Ascents, descents, time/energy share by altitude band and by longitudinal block."""
    z = corr.pos[route.nodes, 2]
    x = corr.pos[route.nodes, 0]
    dt = np.diff(route.times)
    de = np.diff(route.energies)
    dz = np.diff(z)
    length = corr.pos[:, 0].max()
    alt_t, alt_e, blk_t, blk_e = {}, {}, np.zeros(n_blocks), np.zeros(n_blocks)
    for k in range(len(dt)):
        lab = altitude_label(z[k], z[k + 1])
        alt_t[lab] = alt_t.get(lab, 0.0) + dt[k]
        alt_e[lab] = alt_e.get(lab, 0.0) + de[k]
        xm = 0.5 * (x[k] + x[k + 1])
        blk = min(int(xm / length * n_blocks), n_blocks - 1)
        blk_t[blk] += dt[k]
        blk_e[blk] += de[k]
    tot_t, tot_e = dt.sum(), de.sum()
    return dict(
        ascents=int((dz > 0).sum()), descents=int((dz < 0).sum()),
        alt_time_pct={k: 100 * v / tot_t for k, v in alt_t.items()},
        alt_energy_pct={k: 100 * v / tot_e for k, v in alt_e.items()},
        block_time_pct=list(100 * blk_t / tot_t), block_energy_pct=list(100 * blk_e / tot_e))

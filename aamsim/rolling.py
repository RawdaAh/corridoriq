"""Priority-aware rolling-horizon corridor scheduling (Paper 3).

Port of ``AAM_Rolling_Horizon_Priority_ILP_Corrected.ipynb`` (14 Oct 2026).

Requests (reveal, preferred and latest departure, priority) become known over
time. At each roll the controller builds space-time candidate plans (K shortest
spatial paths x feasible departure slots), solves a two-stage lexicographic
binary program (max priority-weighted service, then min weighted ground delay +
flight time at that service level) and commits plans whose departure falls in
the freeze window. Committed missions keep their full reservation footprint.

Baselines: causal FCFS on the same epochs/candidates/freeze rule, and an
offline full-information restricted benchmark.
"""
from __future__ import annotations

import itertools
import math
import time
from collections import defaultdict
from dataclasses import asdict, dataclass, replace

import networkx as nx
import numpy as np
import pandas as pd
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix, vstack

from .corridor import CorridorSpec, build_corridor, mph_to_fps

PRIORITY = {"normal": 1.0, "high": 4.0, "emergency": 12.0}
PRIORITY_PROBS = [0.78, 0.18, 0.04]


@dataclass(frozen=True)
class RHConfig:
    length_ft: float = 500.0
    width_ft: float = 200.0
    height_ft: float = 100.0
    nx: int = 4
    ny: int = 1
    nz: int = 3
    speed_mph: float = 35.0
    energy_j_per_ft: float = 20.42
    dt_min: float = 0.1
    demand_min: float = 30.0
    max_hold_min: float = 8.0
    lead_min: float = 5.0
    roll_min: float = 1.0
    horizon_min: float = 6.0
    freeze_min: float = 2.0
    k_paths: int = 4
    time_limit_sec: float = 10.0
    offline_time_sec: float = 60.0
    mip_gap: float = 0.0
    age_weight: float = 0.0
    ground_weight: float = 1.0
    travel_weight: float = 0.05
    resource_mode: str = "edge"          # edge | edge_node | edge_node_cross
    crossing_tolerance_ft: float = 0.01
    vertical_ground_only: bool = False
    n_requests: int = 50
    pattern: str = "peak"                # uniform | peak | two_peak
    forward_share: float = 0.5
    n_operators: int = 3
    seed: int = 42
    # extension: static wind snapshot changing link travel slots (0 = no wind)
    wind_speed_mph: float = 0.0
    wind_direction_deg: float = 0.0

    def slots(self, minutes: float) -> int:
        value = minutes / self.dt_min
        if not np.isclose(value, round(value), atol=1e-8, rtol=0):
            raise ValueError(f"{minutes} min must be an integer multiple of dt={self.dt_min}")
        return int(round(value))

    def validate(self) -> "RHConfig":
        if not 0 < self.roll_min <= self.freeze_min <= self.horizon_min:
            raise ValueError("Require 0 < roll <= freeze <= planning horizon")
        for key in ["demand_min", "max_hold_min", "lead_min", "roll_min", "horizon_min",
                    "freeze_min"]:
            self.slots(getattr(self, key))
        if self.pattern not in ["uniform", "peak", "two_peak"]:
            raise ValueError("Unknown demand pattern")
        if self.resource_mode not in ["edge", "edge_node", "edge_node_cross"]:
            raise ValueError("Unknown resource mode")
        return self


def segment_distance(a, b, c, d):
    u, v = b - a, d - c
    cands = []
    st = np.linalg.lstsq(np.column_stack([u, -v]), c - a, rcond=None)[0]
    if np.all(st >= 0) and np.all(st <= 1):
        cands.append(np.linalg.norm(a + st[0] * u - c - st[1] * v))
    for p, x, y in [(a, c, d), (b, c, d), (c, a, b), (d, a, b)]:
        w = y - x
        z = np.clip(np.dot(p - x, w) / np.dot(w, w), 0, 1)
        cands.append(np.linalg.norm(p - x - z * w))
    return min(cands)


class RHNetwork:
    """Corridor + rounded travel slots + K-shortest-path cache + resource tokens."""

    def __init__(self, cfg: RHConfig):
        self.cfg = cfg.validate()
        self.corr = build_corridor(CorridorSpec(
            n_long=cfg.nx, n_lat=cfg.ny, n_alt=cfg.nz, module_length=cfg.length_ft,
            module_width=cfg.width_ft, module_height=cfg.height_ft,
            vertical_ground_only=cfg.vertical_ground_only))
        c = self.corr
        self.pos = {i: c.pos[i] for i in range(c.n_nodes)}
        fps = mph_to_fps(cfg.speed_mph)
        lam = np.ones(c.n_arcs)
        if cfg.wind_speed_mph > 0:
            from .wind import from_speed_function
            w = from_speed_function(c, lambda x, y, z, t: np.full_like(x, cfg.wind_speed_mph),
                                    cfg.wind_direction_deg)
            lam = w.link_factors(c, cfg.speed_mph, "continuous")[0][0]
        self.nominal_min = c.arc_len / fps / 60.0 * lam
        self.travel_slots = np.maximum(1, np.ceil(self.nominal_min / cfg.dt_min - 1e-12)).astype(int)
        self.graph = nx.DiGraph()
        self.lookup = {}
        for a, (u, v) in enumerate(zip(c.arc_from, c.arc_to)):
            if np.isfinite(lam[a]):
                self.graph.add_edge(int(u), int(v), weight=int(self.travel_slots[a]))
                self.lookup[(int(u), int(v))] = a
        self.origins, self.destinations = list(c.origins), list(c.destinations)
        self.ground = {i for i in range(c.n_nodes) if c.pos[i, 2] == 0}
        self.cache = {}
        self.crossings = defaultdict(list)
        if cfg.resource_mode == "edge_node_cross":
            edges = [tuple(e) for e in c.edge_nodes.tolist()]
            for e, f in itertools.combinations(edges, 2):
                if set(e) & set(f):
                    continue
                if segment_distance(c.pos[e[0]], c.pos[e[1]], c.pos[f[0]], c.pos[f[1]]) \
                        <= cfg.crossing_tolerance_ft:
                    token = ("cross",) + tuple(sorted((e, f)))
                    self.crossings[e].append(token)
                    self.crossings[f].append(token)

    def paths(self, o, d, k=None):
        k = self.cfg.k_paths if k is None else k
        key = (int(o), int(d), int(k))
        if key not in self.cache:
            g = self.graph.subgraph(set(self.graph) - (self.ground - {int(o), int(d)}))
            try:
                self.cache[key] = [tuple(p) for p in itertools.islice(
                    nx.shortest_simple_paths(g, int(o), int(d), weight="weight"), k)]
            except (nx.NetworkXNoPath, nx.NodeNotFound):
                self.cache[key] = []
        return self.cache[key]

    def trajectory(self, path, departure):
        cfg = self.cfg
        t = int(departure)
        resources, legs = set(), []
        nominal = energy = distance = 0.0
        if cfg.resource_mode != "edge":
            resources.add((("node", path[0]), t))
        for u, v in zip(path[:-1], path[1:]):
            a = self.lookup[(u, v)]
            end = t + int(self.travel_slots[a])
            edge = tuple(sorted((u, v)))
            for slot in range(t, end):
                resources.add((("edge",) + edge, slot))
                for token in self.crossings[edge]:
                    resources.add((token, slot))
            if cfg.resource_mode != "edge":
                resources.add((("node", v), end))
            legs.append((u, v, t, end))
            t = end
            nominal += self.nominal_min[a]
            distance += self.corr.arc_len[a]
            energy += self.corr.arc_len[a] * cfg.energy_j_per_ft
        return dict(node_path=tuple(path), departure_slot=int(departure), arrival_slot=t,
                    travel_slots=t - int(departure), nominal_travel_min=nominal,
                    scheduled_travel_min=(t - int(departure)) * cfg.dt_min,
                    energy_j=energy, distance_ft=distance, resources=frozenset(resources),
                    legs=legs)


REQUEST_COLUMNS = ["request_id", "operator_id", "origin", "destination", "reveal",
                   "preferred", "latest", "priority"]


def rh_demand(cfg: RHConfig, net: RHNetwork) -> pd.DataFrame:
    rng = np.random.default_rng(cfg.seed)
    n = cfg.n_requests
    horizon = cfg.slots(cfg.demand_min)
    if cfg.pattern == "uniform":
        pref = rng.uniform(0, cfg.demand_min, n)
    elif cfg.pattern == "peak":
        pref = rng.normal(0.55 * cfg.demand_min, 0.12 * cfg.demand_min, n)
    else:
        centers = np.where(rng.random(n) < 0.5, 0.30, 0.72) * cfg.demand_min
        pref = rng.normal(centers, 0.08 * cfg.demand_min)
    pref = np.clip(np.rint(pref / cfg.dt_min).astype(int), 0, horizon)
    rows = []
    for rid, p in enumerate(pref, 1):
        forward = rng.random() < cfg.forward_share
        o = int(rng.choice(net.origins if forward else net.destinations))
        d = int(rng.choice(net.destinations if forward else net.origins))
        rows.append(dict(request_id=rid, operator_id=int(rng.integers(1, cfg.n_operators + 1)),
                         origin=o, destination=d,
                         reveal=max(0, int(p) - cfg.slots(cfg.lead_min)), preferred=int(p),
                         latest=int(p) + cfg.slots(cfg.max_hold_min),
                         priority=str(rng.choice(list(PRIORITY), p=PRIORITY_PROBS))))
    return pd.DataFrame(rows, columns=REQUEST_COLUMNS)


def candidates(req, net, now, horizon, reserved):
    plans = []
    for r in req.itertuples():
        for dep in range(max(now, int(r.preferred)), min(int(r.latest), now + horizon - 1) + 1):
            for rank, path in enumerate(net.paths(r.origin, r.destination), 1):
                p = net.trajectory(path, dep)
                if p["resources"] & reserved:
                    continue
                p.update(request_id=int(r.request_id), path_rank=rank,
                         delay_min=(dep - r.preferred) * net.cfg.dt_min,
                         weight=PRIORITY[r.priority]
                         + net.cfg.age_weight * max(0, now - r.reveal) * net.cfg.dt_min)
                plans.append(p)
    return plans


def _master(plans):
    incidence = defaultdict(list)
    for j, p in enumerate(plans):
        incidence[("request", p["request_id"])].append(j)
        for token in p["resources"]:
            incidence[("resource", token)].append(j)
    ii, jj = [], []
    for i, cols in enumerate(incidence.values()):
        ii.extend([i] * len(cols))
        jj.extend(cols)
    return coo_matrix((np.ones(len(ii)), (ii, jj)), shape=(len(incidence), len(plans))).tocsr()


def _valid_binary(x, A, lo, hi):
    if x is None:
        return None
    x = np.asarray(x)
    if x.shape != (A.shape[1],) or not np.isfinite(x).all():
        return None
    z = np.rint(x)
    if np.max(np.abs(x - z), initial=0) > 1e-6 or np.any(z < 0) or np.any(z > 1):
        return None
    az = A @ z
    if np.any(az < lo - 1e-6) or np.any(az > hi + 1e-6):
        return None
    return z.astype(int)


def solve_plans(plans, cfg: RHConfig, time_limit=None):
    """Two-stage lexicographic selection with a feasible greedy fallback."""
    start = time.perf_counter()
    n = len(plans)
    if not n:
        return [], dict(fallback=False, service=0.0, service_upper=0.0, stage1_certified=True,
                        lexicographic_certified=True, solver_sec=0.0, candidates=0)
    A = _master(plans)
    lo, hi = np.full(A.shape[0], -np.inf), np.ones(A.shape[0])
    w = np.array([p["weight"] for p in plans])
    c = np.array([cfg.ground_weight * p["delay_min"] + cfg.travel_weight * p["scheduled_travel_min"]
                  for p in plans])
    greedy = np.zeros(n, dtype=int)
    used, assigned = set(), set()
    for j in sorted(range(n), key=lambda j: (-w[j], c[j], plans[j]["request_id"], j)):
        p = plans[j]
        if p["request_id"] not in assigned and not (p["resources"] & used):
            greedy[j] = 1
            assigned.add(p["request_id"])
            used.update(p["resources"])
    opts = {"time_limit": float(time_limit or cfg.time_limit_sec), "mip_rel_gap": cfg.mip_gap}
    r1 = milp(-w, integrality=np.ones(n), bounds=Bounds(0, 1),
              constraints=LinearConstraint(A, lo, hi), options=opts)
    if r1.status in [2, 3, 4]:
        raise RuntimeError(f"Unexpected stage-1 solver failure: {r1.message}")
    x = _valid_binary(r1.x, A, lo, hi)
    fallback = x is None
    if x is None or (w @ greedy > w @ x + 1e-6):
        x, fallback = greedy, True
    service = float(w @ x)
    A2 = vstack([A, coo_matrix(w.reshape(1, -1))]).tocsr()
    lo2, hi2 = np.r_[lo, service], np.r_[hi, service]
    r2 = milp(c, integrality=np.ones(n), bounds=Bounds(0, 1),
              constraints=LinearConstraint(A2, lo2, hi2), options=opts)
    if r2.status in [2, 3, 4]:
        raise RuntimeError(f"Unexpected stage-2 failure: {r2.message}")
    y = _valid_binary(r2.x, A2, lo2, hi2)
    if y is not None and c @ y <= c @ x + 1e-6:
        x = y
    dual = getattr(r1, "mip_dual_bound", None)
    upper = float(-dual) if dual is not None and np.isfinite(dual) else None
    cert = upper is not None and upper - service <= 1e-5
    dual2 = getattr(r2, "mip_dual_bound", None)
    cert2 = dual2 is not None and np.isfinite(dual2) and c @ x - float(dual2) <= 1e-5
    return [p.copy() for j, p in enumerate(plans) if x[j]], dict(
        fallback=fallback, service=float(w @ x), service_upper=upper,
        stage1_certified=bool(cert), lexicographic_certified=bool(cert and cert2),
        solver_sec=time.perf_counter() - start, candidates=n)


def _signature(p):
    return (p["departure_slot"], p["node_path"])


def rolling(req, net, policy="ilp"):
    """Online controller. ``policy`` is 'ilp' or 'fcfs'."""
    cfg = net.cfg
    statuses = {int(i): "NOT_REVEALED" for i in req.request_id}
    plans, reserved, previous, logs, changes = {}, set(), {}, [], defaultdict(int)
    step, freeze, horizon = cfg.slots(cfg.roll_min), cfg.slots(cfg.freeze_min), cfg.slots(cfg.horizon_min)
    last = int(req.latest.max()) if len(req) else 0
    for now in range(0, (last // step + 1) * step + 1, step):
        wall = time.perf_counter()
        for r in req.itertuples():
            if statuses[r.request_id] == "NOT_REVEALED" and r.reveal <= now:
                statuses[r.request_id] = "QUEUED"
            if statuses[r.request_id] not in ["COMMITTED", "EXPIRED"] and r.latest < now:
                statuses[r.request_id] = "EXPIRED"
        active = req[req.request_id.map(statuses).isin(["QUEUED", "TENTATIVE"])]
        pool = candidates(active, net, now, horizon, reserved)
        if policy == "ilp":
            selected, info = solve_plans(pool, cfg)
        elif policy == "fcfs":
            selected, busy, by_req = [], set(), defaultdict(list)
            for p in pool:
                by_req[p["request_id"]].append(p)
            for r in active.sort_values(["reveal", "preferred", "request_id"]).itertuples():
                for p in sorted(by_req[r.request_id], key=lambda p: (
                        p["departure_slot"], p["scheduled_travel_min"], p["path_rank"])):
                    if not p["resources"] & busy:
                        selected.append(p)
                        busy.update(p["resources"])
                        break
            info = dict(fallback=False, stage1_certified=False, lexicographic_certified=False,
                        solver_sec=0.0, candidates=len(pool),
                        service=sum(p["weight"] for p in selected), service_upper=None)
        else:
            raise ValueError("Unknown policy")
        for rid in active.request_id:
            statuses[int(rid)] = "QUEUED"
        current = {p["request_id"]: _signature(p) for p in selected}
        for rid, sig in previous.items():
            if current.get(rid) != sig:
                changes[rid] += 1
        tentative, commits = {}, 0
        for p in selected:
            rid = p["request_id"]
            if p["departure_slot"] < now + freeze:
                p["commit_slot"] = now
                plans[rid] = p
                reserved.update(p["resources"])
                statuses[rid] = "COMMITTED"
                commits += 1
            else:
                statuses[rid] = "TENTATIVE"
                tentative[rid] = _signature(p)
        previous = tentative
        logs.append(dict(epoch=now, time_min=now * cfg.dt_min, active=len(active),
                         committed=len(plans), new_commits=commits,
                         expired=sum(v == "EXPIRED" for v in statuses.values()),
                         tentative=dict(tentative),
                         wall_sec=time.perf_counter() - wall, **info))
    if any(v not in ["COMMITTED", "EXPIRED"] for v in statuses.values()):
        raise AssertionError("Unresolved requests after final epoch")
    return dict(plans=plans, statuses=statuses, logs=logs, replans=dict(changes), policy=policy)


def offline(req, net):
    horizon = int(req.latest.max()) + 1 if len(req) else 1
    start = time.perf_counter()
    pool = candidates(req, net, 0, horizon, set())
    selected, info = solve_plans(pool, net.cfg, net.cfg.offline_time_sec)
    plans = {p["request_id"]: p for p in selected}
    info["wall_sec"] = time.perf_counter() - start
    return dict(plans=plans, statuses={int(i): "COMMITTED" if i in plans else "EXPIRED"
                                       for i in req.request_id},
                logs=[info], replans={}, policy="offline")


def audit(result, req, net):
    """Rebuild every committed mission from geometry and verify exclusivity/timing."""
    lookup = req.set_index("request_id")
    owner = {}
    cfg = net.cfg
    for rid, p in result["plans"].items():
        r = lookup.loc[rid]
        path, dep = p["node_path"], p["departure_slot"]
        assert path[0] == r.origin and path[-1] == r.destination
        assert len(path) == len(set(path)) and len(path) >= 2
        assert not (set(path[1:-1]) & net.ground)
        assert r.preferred <= dep <= r.latest and dep >= r.reveal
        if result["policy"] != "offline":
            assert r.reveal <= p["commit_slot"] <= dep
            assert dep < p["commit_slot"] + cfg.slots(cfg.freeze_min)
        t, tokens = dep, set()
        if cfg.resource_mode != "edge":
            tokens.add((("node", path[0]), t))
        for u, v in zip(path[:-1], path[1:]):
            a = net.lookup[(u, v)]
            end = t + int(net.travel_slots[a])
            edge = tuple(sorted((u, v)))
            for s in range(t, end):
                tokens.add((("edge",) + edge, s))
                for cross in net.crossings[edge]:
                    tokens.add((cross, s))
            if cfg.resource_mode != "edge":
                tokens.add((("node", v), end))
            t = end
        assert t == p["arrival_slot"]
        assert tokens == set(p["resources"])
        for token in tokens:
            if token in owner:
                raise AssertionError(f"Conflict: {rid}, {owner[token]}, {token}")
            owner[token] = rid
    assert set(result["statuses"]) == set(req.request_id)
    assert all(v in ["COMMITTED", "EXPIRED"] for v in result["statuses"].values())
    return {"passed": True, "requests": len(req), "served": len(result["plans"]),
            "checked_tokens": len(owner)}


def rh_metrics(result, req, net) -> dict:
    cfg = net.cfg
    rows = []
    lookup = req.set_index("request_id")
    for rid, p in result["plans"].items():
        r = lookup.loc[rid]
        rows.append(dict(priority=r.priority, operator_id=int(r.operator_id),
                         delay=p["delay_min"], travel=p["scheduled_travel_min"],
                         nominal=p["nominal_travel_min"], energy=p["energy_j"]))
    s = pd.DataFrame(rows, columns=["priority", "operator_id", "delay", "travel", "nominal", "energy"])
    n = len(req)
    weight = float(sum(PRIORITY[x] for x in s.priority))
    total = float(sum(PRIORITY[x] for x in req.priority))
    logs = result["logs"]
    row = dict(method=result["policy"], requests=n, served=len(s),
               service_rate=len(s) / n if n else np.nan, weighted_service=weight,
               weighted_service_rate=weight / total if total else np.nan,
               mean_delay=s.delay.mean(), p95_delay=s.delay.quantile(0.95) if len(s) else np.nan,
               max_delay=s.delay.max(), scheduled_travel=s.travel.mean(),
               rounding_overhead=(s.travel - s.nominal).mean(),
               total_energy_j=s.energy.sum(), replans=sum(result["replans"].values()),
               mean_wall_sec=float(np.mean([x["wall_sec"] for x in logs])),
               max_wall_sec=float(max(x["wall_sec"] for x in logs)),
               fallback_epochs=int(sum(x.get("fallback", False) for x in logs)),
               stage1_certified_fraction=float(np.mean([x.get("stage1_certified", False)
                                                        for x in logs])))
    for cls in PRIORITY:
        den = int((req.priority == cls).sum())
        row[cls + "_service_rate"] = (s.priority == cls).sum() / den if den else np.nan
    rates = [((s.operator_id == op).sum() / (req.operator_id == op).sum())
             for op in sorted(req.operator_id.unique())]
    row["operator_rate_range"] = max(rates) - min(rates) if rates else np.nan
    return row


def run_case(cfg: RHConfig, methods=("fcfs", "ilp", "offline")):
    """Run all methods on identical demand, audit each, return (table, results, req, net)."""
    net = RHNetwork(cfg)
    req = rh_demand(cfg, net)
    rows, results = [], {}
    for m in methods:
        t0 = time.perf_counter()
        res = offline(req, net) if m == "offline" else rolling(req, net, m)
        audit(res, req, net)
        row = rh_metrics(res, req, net)
        row["total_wall_sec"] = time.perf_counter() - t0
        rows.append(row)
        results[m] = res
    table = pd.DataFrame(rows)
    if "offline" in results and cfg.age_weight == 0:
        up = results["offline"]["logs"][0]["service_upper"]
        table["offline_service_upper"] = up
        table["gap_to_service_upper"] = (up - table.weighted_service) if up is not None else np.nan
    return table, results, req, net


def experiment_design(base: RHConfig, quick: bool = False):
    """One-factor-at-a-time grid of the notebook (subset when quick=True)."""
    cases = [("base", base)]
    loads = [100, 200, 300] if quick else [100, 150, 200, 250, 300, 350, 400]
    cases += [(f"load_{n}", replace(base, n_requests=n)) for n in loads]
    cases += [(f"demand_{p}", replace(base, pattern=p)) for p in ["uniform", "two_peak"]]
    cases += [(f"freeze_{f}", replace(base, freeze_min=f)) for f in [1.0, 4.0]]
    cases += [(f"horizon_{h}", replace(base, horizon_min=h)) for h in [4.0, 10.0]]
    cases += [(f"paths_{k}", replace(base, k_paths=k)) for k in [1, 8]]
    if not quick:
        cases += [(f"roll_{r}", replace(base, roll_min=r)) for r in [0.5, 2.0]]
        cases += [(f"dt_{d}", replace(base, dt_min=d)) for d in [0.05, 0.25, 1.0]]
        cases += [(m, replace(base, resource_mode=m)) for m in ["edge_node", "edge_node_cross"]]
        cases += [("vertical_ground", replace(base, vertical_ground_only=True))]
        cases += [(f"aging_{a}", replace(base, age_weight=a)) for a in [0.12, 0.5]]
        cases += [(f"lead_{l}", replace(base, lead_min=l)) for l in [0.0, 10.0]]
        cases += [(f"direction_{f}", replace(base, forward_share=f)) for f in [0.0, 1.0]]
        cases += [("larger_network", replace(base, nx=8, ny=2, nz=3))]
    cases += [(f"wind_{d:g}", replace(base, wind_speed_mph=15.0, wind_direction_deg=d))
              for d in [0.0, 180.0]]
    for _, c in cases:
        c.validate()
    return cases


def config_dict(cfg: RHConfig) -> dict:
    return asdict(cfg)

"""Shared-Intent Reservation Engine (SIRE) - TRB paper (Paper 2).

Sequential pre-departure registry. Requests are processed by intended dispatch
slot (ties broken deterministically by operator id, then trip id). When a request
reaches its turn the operator reads the published occupancy ledger, computes its
own minimum-travel-time route on the time-expanded network over the open
link-time resources (Eqs. 10-12) and submits it; the registry commits it
atomically (Eq. 13) or rejects the request.

The route search is a label-setting search on the time-expanded state space
(node, slot). Because every transition moves forward in time by at least one
slot, it is evaluated layer by layer with numpy, which is equivalent to the
Dijkstra search of the original notebook but much faster.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .corridor import Corridor
from .energy import VehicleModel, baseline_link_costs

INF = np.inf

try:  # optional compiled kernel; the pure-numpy search is used otherwise
    import numba as _nb
    _HAVE_NUMBA = True
except ImportError:  # pragma: no cover
    _HAVE_NUMBA = False


def _route_kernel_py(ptr, order, af, at, key, tau, q, blocked, origin, dest, t0, H, noland):
    """Label-setting search on (node, slot) layers. Compiled with numba if present."""
    N = ptr.shape[0] - 1
    best = np.full((H + 1, N), np.inf)
    pred = np.full((H + 1, N), -1, np.int64)
    pred_t = np.full((H + 1, N), -1, np.int64)
    best[t0, origin] = 0.0
    best_dest = np.inf
    best_t = -1
    frontier_end = t0
    for t in range(t0, H):
        if t > frontier_end:
            break
        c = best[t, dest]
        if c < best_dest:
            best_dest = c
            best_t = t
        for u in range(N):
            cu = best[t, u]
            if cu == np.inf or cu >= best_dest or u == dest:
                continue
            for p in range(ptr[u], ptr[u + 1]):
                a = order[p]
                if noland[at[a]] and at[a] != dest:
                    continue
                qa = q[t, a]
                ta = t + qa
                if ta >= H:
                    continue
                k = key[a]
                free = True
                for s in range(t, ta):
                    if blocked[k, s]:
                        free = False
                        break
                if not free:
                    continue
                cand = cu + tau[t, a]
                v = at[a]
                if cand < best[ta, v]:
                    best[ta, v] = cand
                    pred[ta, v] = a
                    pred_t[ta, v] = t
                    if ta > frontier_end:
                        frontier_end = ta
    for t in range(t0, H + 1):
        if best[t, dest] < best_dest:
            best_dest = best[t, dest]
            best_t = t
    return best_dest, best_t, pred, pred_t


_route_kernel = _nb.njit(cache=True)(_route_kernel_py) if _HAVE_NUMBA else _route_kernel_py


@dataclass
class SireConfig:
    slot_sec: float = 60.0            # reservation slot duration (Delta t)
    dispatch_window: int = 60         # slots in which trips are dispatched
    route_horizon: int = 120          # latest allowable arrival slot (exclusive)
    conflict_rule: str = "undirected"  # 'undirected' (shared) | 'directed' (separated)
    tie_break: str = "operator"       # 'operator' | 'random' | 'rotating'
    dispatch_flex: int = 0            # extra slots a trip may wait on the ground
    ground_transit: bool = True       # False: no landing at intermediate ground nodes
    seed: int = 0


@dataclass
class Trip:
    trip_id: int
    operator: int
    origin: int
    destination: int
    dispatch: int


@dataclass
class SireResult:
    trips: pd.DataFrame                 # one row per request with outcome
    legs: list = field(default_factory=list)  # per accepted trip: list of (arc, t_dep, t_arr)
    config: SireConfig | None = None

    def kpis(self, window: int = 5) -> dict:
        return sire_kpis(self.trips, window)


# ----------------------------------------------------------------------------- demand

def generate_trips(corr: Corridor, n_trips: int, n_operators: int = 3,
                   window: int = 60, pattern: str = "uniform", forward_share: float = 1.0,
                   seed: int = 0) -> list[Trip]:
    """Balanced demand generator of the TRB paper.

    Demand per operator is split as evenly as possible; within each operator the
    forward share goes entry->exit face, the rest exit->entry. Each direction is
    spread evenly over its OD pairs with any remainder assigned to the first pairs
    in index order. Dispatch slots are uniform on [0, window-1] or, for
    ``pattern='peak'``, normal around window/2 with sd window/4 (rounded, clipped).
    """
    rng = np.random.default_rng(seed)
    fwd_pairs = [(o, d) for o in corr.origins for d in corr.destinations]
    rev_pairs = [(d, o) for o in corr.origins for d in corr.destinations]
    per_op = [n_trips // n_operators + (1 if m < n_trips % n_operators else 0)
              for m in range(n_operators)]
    trips = []
    tid = 0
    for m, dm in enumerate(per_op, start=1):
        n_fwd = int(round(forward_share * dm))
        ods = []
        for pairs, n in ((fwd_pairs, n_fwd), (rev_pairs, dm - n_fwd)):
            base, rem = divmod(n, len(pairs))
            for idx, od in enumerate(pairs):
                ods.extend([od] * (base + (1 if idx < rem else 0)))
        if pattern == "uniform":
            times = rng.integers(0, window, size=len(ods))
        elif pattern == "peak":
            times = rng.normal(window // 2, max(1.0, window * 0.25), size=len(ods))
            times = np.clip(np.rint(times), 0, window - 1).astype(int)
        else:
            raise ValueError("pattern must be 'uniform' or 'peak'")
        rng.shuffle(times)
        for (o, d), t in zip(ods, times):
            tid += 1
            trips.append(Trip(tid, m, int(o), int(d), int(t)))
    return trips


# ------------------------------------------------------------------------ the engine

class SireEngine:
    """Reservation ledger + operator route search on one corridor."""

    def __init__(self, corr: Corridor, veh: VehicleModel | None = None,
                 cfg: SireConfig | None = None, link_time_sec=None):
        self.corr = corr
        self.veh = veh or VehicleModel()
        self.cfg = cfg or SireConfig()
        H = self.cfg.route_horizon
        if link_time_sec is None:
            t, _ = baseline_link_costs(corr, self.veh)
            link_time_sec = np.broadcast_to(t, (H + 1, corr.n_arcs))
        lt = np.asarray(link_time_sec, dtype=float)
        if lt.ndim == 1:
            lt = np.broadcast_to(lt, (H + 1, corr.n_arcs))
        self.tau = lt                                       # (H+1, A) seconds by departure slot
        self.q = np.maximum(1, np.ceil(lt / self.cfg.slot_sec - 1e-9)).astype(int)  # Eq. 7
        if self.cfg.conflict_rule == "undirected":
            self.key = corr.arc_edge                        # Eq. 6
            n_keys = corr.n_edges
        elif self.cfg.conflict_rule == "directed":
            self.key = np.arange(corr.n_arcs)
            n_keys = corr.n_arcs
        else:
            raise ValueError("conflict_rule must be 'undirected' or 'directed'")
        self.qmax = int(self.q.max())
        self.blocked = np.zeros((n_keys, H + self.qmax + 2), dtype=bool)
        self.H = H
        self.use_numba = True

    # ---------------------------------------------------------------- route search
    def route(self, origin: int, destination: int, t0: int):
        """Minimum-travel-time route departing exactly at slot t0 (no airborne wait).

        Returns (cost_sec, legs) with legs = [(arc, t_dep, t_arr), ...] or None.
        """
        if not 0 <= t0 < self.H:
            return None
        if _HAVE_NUMBA and self.use_numba:
            if not hasattr(self, "_csr"):
                order = np.argsort(self.corr.arc_from, kind="stable")
                ptr = np.zeros(self.corr.n_nodes + 1, dtype=np.int64)
                np.add.at(ptr, self.corr.arc_from + 1, 1)
                self._csr = (np.cumsum(ptr), order.astype(np.int64))
                self._noland = (self.corr.ground & (not self.cfg.ground_transit)).astype(np.bool_)
                self._tau_c = np.ascontiguousarray(self.tau)
                self._q_c = np.ascontiguousarray(self.q).astype(np.int64)
            ptr, order = self._csr
            cost, tdest, pred, pred_t = _route_kernel(
                ptr, order, self.corr.arc_from, self.corr.arc_to, self.key.astype(np.int64),
                self._tau_c, self._q_c, self.blocked, int(origin), int(destination), int(t0), self.H,
                self._noland)
            if tdest < 0:
                return None
            return float(cost), self._unwind(origin, destination, t0, tdest, pred, pred_t)
        return self._route_numpy(origin, destination, t0)

    def _unwind(self, origin, destination, t0, tdest, pred, pred_t):
        legs = []
        node, t = destination, tdest
        while not (node == origin and t == t0):
            a, s = int(pred[t, node]), int(pred_t[t, node])
            legs.append((a, s, t))
            node, t = int(self.corr.arc_from[a]), s
        legs.reverse()
        return legs

    def _route_numpy(self, origin: int, destination: int, t0: int):
        corr, H = self.corr, self.H
        N = corr.n_nodes
        best = np.full((H + 1, N), INF)
        pred = np.full((H + 1, N), -1, dtype=np.int64)
        pred_t = np.full((H + 1, N), -1, dtype=np.int64)
        best[t0, origin] = 0.0
        best_dest, best_dest_t = INF, -1
        frontier_end = t0
        af, at, key = corr.arc_from, corr.arc_to, self.key
        for t in range(t0, H):
            if t > frontier_end:
                break
            row = best[t]
            if np.isfinite(row[destination]) and row[destination] < best_dest:
                best_dest, best_dest_t = row[destination], t
            act = np.isfinite(row) & (row < best_dest)
            act[destination] = False
            if not act.any():
                continue
            arcs = np.flatnonzero(act[af])
            if not self.cfg.ground_transit:
                land = corr.ground[at[arcs]] & (at[arcs] != destination)
                arcs = arcs[~land]
            if arcs.size == 0:
                continue
            q = self.q[t, arcs]
            if self.qmax == 1:
                ok = ~self.blocked[key[arcs], t]
            else:
                ok = np.ones(arcs.size, dtype=bool)
                for s in range(self.qmax):
                    m = q > s
                    ok[m] &= ~self.blocked[key[arcs[m]], t + s]
            arcs, q = arcs[ok], q[ok]
            t_arr = t + q
            keep = t_arr < H
            arcs, t_arr = arcs[keep], t_arr[keep]
            if arcs.size == 0:
                continue
            cand = row[af[arcs]] + self.tau[t, arcs]
            # scatter-min per (arrival slot, node), remembering the winning arc
            tgt = t_arr * N + at[arcs]
            order = np.lexsort((cand, tgt))
            tgt_s, first = np.unique(tgt[order], return_index=True)
            win = order[first]
            flat_best = best.reshape(-1)
            better = cand[win] < flat_best[tgt_s]
            upd = tgt_s[better]
            flat_best[upd] = cand[win][better]
            pred.reshape(-1)[upd] = arcs[win][better]
            pred_t.reshape(-1)[upd] = t
            if upd.size:
                frontier_end = max(frontier_end, int(t_arr[win][better].max()))
        # final check of the last layers
        col = best[:, destination]
        if np.isfinite(col).any():
            tt = int(np.argmin(col))  # earliest slot among equal minimum cost
            if col[tt] < best_dest:
                best_dest, best_dest_t = col[tt], tt
        if not np.isfinite(best_dest):
            return None
        return float(best_dest), self._unwind(origin, destination, t0, best_dest_t, pred, pred_t)

    def is_free(self, legs) -> bool:
        for a, s, e in legs:
            if self.blocked[self.key[a], s:e].any():
                return False
        return True

    def commit(self, legs):
        for a, s, e in legs:
            self.blocked[self.key[a], s:e] = True

    # ------------------------------------------------------------- sequential process
    def run(self, trips: list[Trip]) -> SireResult:
        cfg = self.cfg
        rng = np.random.default_rng(cfg.seed + 104729)
        df = pd.DataFrame([t.__dict__ for t in trips])
        if df.empty:
            return SireResult(df.assign(accepted=[], travel_sec=[], arrival=[], order=[]), [], cfg)
        if cfg.tie_break == "operator":
            df = df.sort_values(["dispatch", "operator", "trip_id"])
        elif cfg.tie_break == "random":
            df = df.assign(_r=rng.random(len(df))).sort_values(["dispatch", "_r"]).drop(columns="_r")
        elif cfg.tie_break == "rotating":
            n_op = df.operator.max()
            df = df.assign(_r=(df.operator - 1 - df.dispatch) % n_op)
            df = df.sort_values(["dispatch", "_r", "trip_id"]).drop(columns="_r")
        else:
            raise ValueError("unknown tie_break")
        accepted, travel, arrival, dep, legs_all = [], [], [], [], []
        for r in df.itertuples():
            res = None
            for t0 in range(r.dispatch, r.dispatch + cfg.dispatch_flex + 1):
                res = self.route(r.origin, r.destination, t0)
                if res is not None:
                    break
            if res is not None and self.is_free(res[1]):   # registry validation
                self.commit(res[1])
                accepted.append(True)
                travel.append(res[0])
                arrival.append(res[1][-1][2])
                dep.append(res[1][0][1])
                legs_all.append(res[1])
            else:
                accepted.append(False)
                travel.append(np.nan)
                arrival.append(-1)
                dep.append(-1)
                legs_all.append(None)
        df = df.assign(order=np.arange(len(df)), accepted=accepted, travel_sec=travel,
                       departure=dep, arrival=arrival).reset_index(drop=True)
        return SireResult(df, legs_all, cfg)


def run_sire(corr: Corridor, trips, cfg: SireConfig | None = None,
             veh: VehicleModel | None = None, link_time_sec=None) -> SireResult:
    return SireEngine(corr, veh, cfg, link_time_sec).run(trips)


# ---------------------------------------------------------------------------- metrics

def gini(values) -> float:
    """Gini coefficient of operator success rates (TRB Eq. 18)."""
    a = np.asarray(values, dtype=float)
    if a.size == 0 or a.mean() == 0:
        return 0.0
    return float(np.abs(a[:, None] - a[None, :]).sum() / (2 * a.size ** 2 * a.mean()))


def cascade_indicator(rejected, window: int = 5) -> float:
    """Share of rejections followed by another rejection within ``window`` (Eq. 19)."""
    f = np.asarray(rejected, dtype=int)
    if f.sum() == 0:
        return 0.0
    n = len(f)
    c = np.concatenate([[0], np.cumsum(f)])
    idx = np.flatnonzero(f)
    nxt = c[np.minimum(idx + 1 + window, n)] - c[idx + 1]
    return float((nxt > 0).sum() / f.sum())


def sire_kpis(trips: pd.DataFrame, window: int = 5) -> dict:
    d = len(trips)
    if d == 0:
        return dict(demand=0, accepted=0, success_rate=np.nan, tt_acc_sec=np.nan,
                    gini=np.nan, cascade=np.nan)
    acc = trips.accepted.to_numpy()
    op_sr = trips.groupby("operator").accepted.mean()
    out = dict(demand=d, accepted=int(acc.sum()), success_rate=float(acc.mean()),
               tt_acc_sec=float(trips.travel_sec[acc].mean()) if acc.any() else np.nan,
               gini=gini(op_sr.values),
               cascade=cascade_indicator(~trips.sort_values("order").accepted.to_numpy(), window))
    for op, v in op_sr.items():
        out[f"sr_op{op}"] = float(v)
    return out

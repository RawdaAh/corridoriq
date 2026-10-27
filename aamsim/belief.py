"""Value-of-wind-information experiment (bridge between Papers 1 and 2).

Mirrors the GeoWind2Plan evaluation protocol: every mission is planned with
four planning winds and every selected route is then flown in the true wind.

    ground_truth  - plan with the true field (perfect information reference)
    forecast      - biased, rotated, lagged and noisy version of the truth
    profile       - height-only profile from one sounding (uniform per layer)
    wind_agnostic - zero wind

Outputs per mission and method: planned and realized time/energy, energy regret
versus ground-truth planning, and the reservation drift - how far the realized
link entry times move away from the link-time slots a SIRE registry would have
reserved from the plan. The drift answers: how large a time buffer must a
reservation protocol hold when operators plan with imperfect wind?
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .corridor import Corridor, CorridorSpec, build_corridor
from .energy import J_PER_KWH, VehicleModel, wind_link_costs
from .pareto import evaluate_route, pareto_routes
from .wind import planning_beliefs, synthetic_urban_wind

METHODS = ["ground_truth", "forecast", "profile", "wind_agnostic"]
METHOD_LABELS = {"ground_truth": "Ground truth", "forecast": "Forecast",
                 "profile": "Profile wind", "wind_agnostic": "Wind-agnostic"}


def belief_corridor(n_long: int = 10, n_lat: int = 2) -> Corridor:
    """Same layers and lanes as the TRR corridor, at 3,000 ft longitudinal resolution
    so that the spatial structure of the wind matters."""
    return build_corridor(CorridorSpec(n_long=n_long, n_lat=n_lat, n_alt=3,
                                       module_length=30000.0 / n_long, module_width=100.0,
                                       z_levels=(0.0, 200.0, 300.0, 400.0)))


@dataclass
class Mission:
    index: int
    origin: int
    destination: int
    depart_bin: int
    wind_dir: float        # direction the wind blows toward (deg, from +x)
    heading: float         # 0 = entry->exit, 180 = exit->entry
    seed: int

    @property
    def wind_class(self) -> str:
        rel = abs((self.wind_dir - self.heading + 180) % 360 - 180)
        return "tailwind" if rel <= 45 else ("headwind" if rel >= 135 else "crosswind")


def sample_missions(corr: Corridor, n: int, seed: int = 0, window_bins: int = 10):
    rng = np.random.default_rng(seed)
    out = []
    for k in range(n):
        fwd = rng.random() < 0.5
        o = int(rng.choice(corr.origins if fwd else corr.destinations))
        d = int(rng.choice(corr.destinations if fwd else corr.origins))
        out.append(Mission(k, o, d, int(rng.integers(0, window_bins)),
                           float(rng.uniform(0, 360)), 0.0 if fwd else 180.0,
                           int(rng.integers(0, 2**31 - 1))))
    return out


def select_route(routes, rule: str = "energy"):
    if not routes:
        return None
    if rule == "energy":
        return min(routes, key=lambda r: (r.energy_j, r.travel_sec))
    if rule == "time":
        return min(routes, key=lambda r: (r.travel_sec, r.energy_j))
    # balanced: minimise normalised sum
    t = np.array([r.travel_sec for r in routes])
    e = np.array([r.energy_j for r in routes])
    score = t / t.min() + e / e.min()
    return routes[int(np.argmin(score))]


def slot_drift(planned_times, realized_times, depart_sec, slot_sec):
    """Max absolute slot shift of link entry times and share of links whose realized
    occupancy leaves the reserved slot."""
    p = np.floor((depart_sec + np.asarray(planned_times[:-1])) / slot_sec)
    r = np.floor((depart_sec + np.asarray(realized_times[:-1])) / slot_sec)
    shift = np.abs(r - p)
    return float(shift.max()) if shift.size else 0.0, float((shift > 0).mean()) if shift.size else 0.0


def run_mission(corr, veh, m: Mission, wind_speed: float, rule: str = "energy",
                model: str = "paper", slot_sec: float = 60.0, bin_sec: float = 180.0,
                n_bins: int = 20, horizon_sec: float = 3600.0):
    truth = synthetic_urban_wind(corr, m.seed, ref_speed_mph=wind_speed,
                                 direction_deg=m.wind_dir, n_bins=n_bins, bin_sec=bin_sec)
    beliefs = planning_beliefs(truth, corr, m.seed)
    truth_costs = wind_link_costs(corr, veh, truth, model)
    depart = m.depart_bin * bin_sec
    rows, routes = [], {}
    for name in METHODS:
        w = beliefs[name]
        R = pareto_routes(corr, veh, w, m.origin, m.destination, depart_sec=depart,
                          horizon_sec=horizon_sec, model=model)
        r = select_route(R, rule)
        if r is None:
            rows.append(dict(method=name, feasible=False))
            continue
        times, energies = evaluate_route(corr, veh, truth, r.arcs, depart, model, truth_costs,
                                         takeoff_j=r.energies[0])
        ok = bool(np.isfinite(times[-1]))
        max_shift, frac_shift = slot_drift(r.times, times, depart, slot_sec) if ok else (np.nan, np.nan)
        z = corr.pos[r.nodes, 2]
        dist_ft = float(corr.arc_len[r.arcs].sum())
        rows.append(dict(
            method=name, feasible=ok, n_pareto=len(R),
            planned_min=r.travel_sec / 60, planned_kwh=r.energy_j / J_PER_KWH,
            realized_min=times[-1] / 60 if ok else np.nan,
            realized_kwh=energies[-1] / J_PER_KWH if ok else np.nan,
            wh_per_km=energies[-1] / 3600 / (dist_ft * 0.0003048) if ok else np.nan,
            arrival_error_sec=(times[-1] - r.travel_sec) if ok else np.nan,
            max_slot_shift=max_shift, share_links_shifted=frac_shift,
            mean_alt_ft=float(np.average(0.5 * (z[1:] + z[:-1]), weights=np.diff(r.times)))
            if len(z) > 1 else np.nan,
            ground_speed_mph=(dist_ft / times[-1]) * 3600 / 5280 if ok else np.nan))
        routes[name] = dict(nodes=[int(n) for n in r.nodes], planned_t=list(map(float, r.times)),
                            realized_t=list(map(float, times)) if ok else None)
    ref = next((x for x in rows if x["method"] == "ground_truth" and x.get("feasible")), None)
    for x in rows:
        x.update(mission=m.index, wind_class=m.wind_class, wind_dir=m.wind_dir,
                 heading=m.heading, depart_bin=m.depart_bin, wind_speed=wind_speed)
        if ref is not None and x.get("feasible"):
            x["extra_energy_pct"] = 100 * (x["realized_kwh"] / ref["realized_kwh"] - 1)
            x["extra_time_pct"] = 100 * (x["realized_min"] / ref["realized_min"] - 1)
    return rows, routes, truth, beliefs

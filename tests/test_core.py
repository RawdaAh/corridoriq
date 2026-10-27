"""Correctness checks: geometry, wind projection, reproduction of published values,
reservation exclusivity and the rolling-horizon audit."""
import numpy as np
import pandas as pd
import pytest

from aamsim.corridor import build_corridor, trb_sire_corridor, trr_wind_corridor
from aamsim.energy import VehicleModel
from aamsim.pareto import pareto_routes
from aamsim.sire import SireConfig, SireEngine, cascade_indicator, generate_trips, gini, run_sire
from aamsim.wind import trr_altitude_wind, zero_wind


def test_module_has_28_links_minus_ground():
    c = build_corridor(n_long=1, n_lat=1, n_alt=1)
    # 8 corners -> 28 pairs, minus the 6 pairs that are both on the ground
    assert c.n_edges == 22
    assert c.n_arcs == 44


def test_vertical_ground_only():
    c = build_corridor(n_long=1, n_lat=1, n_alt=1, vertical_ground_only=True)
    g = c.ground
    for u, v in c.edge_nodes:
        if g[u] or g[v]:
            assert np.allclose(c.pos[u, :2], c.pos[v, :2])


def test_sire_geometry_matches_paper():
    s = trb_sire_corridor().summary()
    assert s["nodes"] == 465 and s["length_ft"] == 30000 and s["ceiling_ft"] == 400


@pytest.mark.parametrize("direction,time_min,energy_kwh", [
    (None, 9.74, 0.1702), (0, 5.50, 0.0539), (180, 16.29, 0.4775)])
def test_trr_table1_reproduction(direction, time_min, energy_kwh):
    c = trr_wind_corridor()
    w = zero_wind(c) if direction is None else trr_altitude_wind(c, direction)
    R = pareto_routes(c, VehicleModel(), w, c.node(0, 0, 0), c.node(5, 0, 0))
    assert abs(np.mean([r.travel_min for r in R]) - time_min) < 0.01
    assert abs(np.mean([r.energy_kwh for r in R]) - energy_kwh) < 0.0005


def test_crosswind_equals_no_wind_under_paper_rule():
    c = trr_wind_corridor()
    R = pareto_routes(c, VehicleModel(), trr_altitude_wind(c, 90), c.node(0, 0, 0), c.node(5, 0, 0))
    assert abs(R[0].travel_min - 9.74) < 0.01


def test_pareto_front_is_nondominated():
    c = trr_wind_corridor()
    R = pareto_routes(c, VehicleModel(), trr_altitude_wind(c, 45), c.node(0, 0, 0), c.node(5, 0, 0))
    for a in R:
        for b in R:
            if a is not b:
                assert not (a.travel_sec <= b.travel_sec and a.energy_j <= b.energy_j)


def test_sire_reservations_never_overlap():
    c = trb_sire_corridor()
    eng = SireEngine(c, cfg=SireConfig())
    res = eng.run(generate_trips(c, 800, seed=3))
    owner = {}
    for legs in res.legs:
        if not legs:
            continue
        for a, s, e in legs:
            for t in range(s, e):
                k = (int(c.arc_edge[a]), t)
                assert k not in owner
                owner[k] = True
    assert 0.7 < res.kpis()["success_rate"] < 0.9


def test_sire_numpy_and_numba_agree():
    c = trb_sire_corridor()
    trips = generate_trips(c, 300, seed=5)
    a = SireEngine(c).run(trips)
    e = SireEngine(c)
    e.use_numba = False
    b = e.run(trips)
    assert (a.trips.accepted == b.trips.accepted).all()
    assert np.allclose(a.trips.travel_sec.fillna(0), b.trips.travel_sec.fillna(0))


def test_sire_published_reference_replicate():
    """TRB notebook R0 check: D = 2000, replicate 1 gives SR = 0.3645."""
    c = trb_sire_corridor()
    r = run_sire(c, generate_trips(c, 2000, seed=1), SireConfig())
    assert abs(r.kpis()["success_rate"] - 0.3645) < 0.01


def test_metrics():
    assert gini([1, 1, 1]) == 0
    assert abs(gini([1, 0, 0]) - 2 / 3) < 1e-9
    assert cascade_indicator([0, 1, 1, 0, 0, 0, 0, 0, 0, 1], window=5) == pytest.approx(1 / 3)


def test_rolling_horizon_tiny_case_is_audited():
    from dataclasses import replace
    from aamsim.rolling import RHConfig, run_case
    cfg = replace(RHConfig(), n_requests=8, nx=2, ny=1, nz=1, demand_min=5.0, time_limit_sec=5,
                  offline_time_sec=10)
    table, res, req, net = run_case(cfg)   # run_case audits every method
    assert set(table.method) == {"fcfs", "ilp", "offline"}
    up = table.offline_service_upper.iloc[0]
    if up is not None and np.isfinite(up):
        assert (table.weighted_service <= up + 1e-6).all()

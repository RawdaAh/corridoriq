"""Real-city experiment: the corridor over downtown Dallas (OpenStreetMap buildings).

Every mission is planned with four planning winds and flown in the building-aware
ground-truth wind, as in GeoWind2Plan. Missions run start->goal and goal->start
under sampled background wind directions; results are reported for tailwind,
headwind and crosswind missions.

    python run_geo_mission.py --missions 64 --speeds 2 4 8
    python run_geo_mission.py --direction 0 --speed 8          # one mission, printed

Data: data/osm/dallas.osm.gz (OpenStreetMap contributors, ODbL). Download another
area with ``--fetch`` and ``--bbox min_lon min_lat max_lon max_lat``.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from aamsim.geo import DALLAS, build_geo_corridor, fetch_osm, load_buildings, nearest_portal
from aamsim.geo_plan import LABELS, METHODS, run_geo_mission
from aamsim.plotting import PALETTE, panel_label, plt, save
from aamsim.urbanwind import WindParams


def setup(osm, start, goal):
    city = load_buildings(osm, name="city")
    gc = build_geo_corridor(city, start, goal)
    return gc


def table(df):
    out = []
    for cls in ["all", "tailwind", "headwind", "crosswind"]:
        sub = df if cls == "all" else df[df.wind_class == cls]
        for m in METHODS:
            s = sub[(sub.method == m) & sub.feasible.astype(bool)]
            out.append(dict(missions=cls, method=LABELS[m], n=len(s), wh_per_km=s.wh_per_km.mean(),
                            sd=s.wh_per_km.std(), extra_energy_pct=s.extra_energy_pct.mean(),
                            saving_vs_agnostic_pct=s.saving_vs_agnostic_pct.mean(),
                            saving_recovered_pct=s.saving_recovered_pct.median(),
                            arrival_error_s=s.arrival_error_s.abs().mean(),
                            mean_alt_m=s.mean_alt_m.mean(), ground_speed_ms=s.ground_speed_ms.mean()))
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--osm", default="data/osm/dallas.osm.gz")
    ap.add_argument("--fetch", action="store_true", help="download --bbox from the OSM API into --osm")
    ap.add_argument("--bbox", type=float, nargs=4, default=DALLAS["bbox"])
    ap.add_argument("--start", type=float, nargs=2, default=DALLAS["start"])
    ap.add_argument("--goal", type=float, nargs=2, default=DALLAS["goal"])
    ap.add_argument("--missions", type=int, default=48, help="wind directions sampled per speed")
    ap.add_argument("--speeds", type=float, nargs="+", default=[2.0, 4.0, 8.0])
    ap.add_argument("--direction", type=float, default=None)
    ap.add_argument("--speed", type=float, default=4.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/geo")
    a = ap.parse_args()
    if a.fetch:
        Path(a.osm).parent.mkdir(parents=True, exist_ok=True)
        fetch_osm(a.bbox, a.osm)
    gc = setup(a.osm, tuple(a.start), tuple(a.goal))
    o, d = nearest_portal(gc, "origin"), nearest_portal(gc, "destination")
    print(f"Corridor {gc.grid['L']:.0f} m, {gc.corridor.n_nodes} nodes, {gc.corridor.n_edges} links "
          f"({gc.removed_edges} links removed by {len(gc.city.buildings)} buildings)")
    if a.direction is not None:
        rows, _, _ = run_geo_mission(gc, WindParams(a.speed, a.direction), o, d)
        print(pd.DataFrame(rows)[["method", "wind_class", "realized_wh", "extra_energy_pct",
                                  "saving_vs_agnostic_pct", "realized_s", "mean_alt_m"]].round(2).to_string())
        return
    rng = np.random.default_rng(a.seed)
    rows = []
    for s in a.speeds:
        for k in range(a.missions):
            direction = float(rng.uniform(0, 360))
            fwd = k % 2 == 0
            r, _, _ = run_geo_mission(gc, WindParams(s, direction), o if fwd else d, d if fwd else o)
            for x in r:
                x["mission"] = k
            rows += r
        print(f"  speed {s:g} m/s: {a.missions} missions done", flush=True)
    df = pd.DataFrame(rows)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "missions.csv", index=False)
    tab = table(df)
    tab.to_csv(out / "summary.csv", index=False)
    pd.set_option("display.width", 180)
    print(tab.round(2).to_string(index=False))
    fig, ax = plt.subplots(figsize=(5.5, 2.8))
    classes = ["tailwind", "headwind", "crosswind"]
    x = np.arange(3)
    for k, m in enumerate(["ground_truth", "forecast", "profile"]):
        vals = [df[(df.wind_class == c) & (df.method == m)].saving_vs_agnostic_pct.mean() for c in classes]
        ax.bar(x + (k - 1) * 0.26, vals, 0.24, color=PALETTE[k], label=LABELS[m])
    ax.set_xticks(x, [c.capitalize() for c in classes])
    ax.set(ylabel="Energy saved vs wind-agnostic (%)")
    ax.grid(axis="x", visible=False)
    ax.legend(fontsize=7)
    save(fig, out / "figures", "geo_energy_saving")
    print(f"Results written to {out}/")


if __name__ == "__main__":
    main()

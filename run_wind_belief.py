"""Bridge experiment - value of wind information (GeoWind2Plan-style protocol).

Each mission is planned with four planning winds (ground truth, forecast,
height-only profile, wind-agnostic). Every selected route is then flown in the
ground-truth wind. Results are reported for all missions and separately for
tailwind, headwind and crosswind missions.

    python run_wind_belief.py --missions 256 --wind-speed 15 --workers 8
    python run_wind_belief.py --mission 0                # one mission, printed
    python run_wind_belief.py --wind-speed 8 --select time

Output: results/wind_belief/missions_<speed>mph_<rule>.csv and figures/.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd

from aamsim.belief import (METHOD_LABELS, METHODS, belief_corridor, run_mission,
                           sample_missions)
from aamsim.energy import VehicleModel
from aamsim.plotting import PALETTE, panel_label, plt, save


def _work(m, speed, rule, model, lanes):
    corr = belief_corridor(n_lat=lanes)
    rows, _, _, _ = run_mission(corr, VehicleModel(), m, speed, rule, model)
    return rows


def summary_table(df):
    """GeoWind2Plan-style table: Wh/km mean +- sd, extra energy vs ground truth, altitude, speed."""
    out = []
    for cls in ["all", "tailwind", "headwind", "crosswind"]:
        sub = df if cls == "all" else df[df.wind_class == cls]
        for m in METHODS:
            s = sub[(sub.method == m) & sub.feasible.astype(bool)]
            out.append(dict(missions=cls, method=METHOD_LABELS[m], n=len(s),
                            wh_per_km=s.wh_per_km.mean(), wh_per_km_sd=s.wh_per_km.std(),
                            extra_energy_pct=s.extra_energy_pct.mean(),
                            extra_time_pct=s.extra_time_pct.mean(),
                            arrival_error_sec=s.arrival_error_sec.abs().mean(),
                            p95_arrival_error_sec=s.arrival_error_sec.abs().quantile(0.95),
                            max_slot_shift=s.max_slot_shift.mean(),
                            share_links_shifted=s.share_links_shifted.mean(),
                            mean_alt_ft=s.mean_alt_ft.mean(),
                            ground_speed_mph=s.ground_speed_mph.mean()))
    return pd.DataFrame(out)


def figures(df, out, tag):
    classes = ["tailwind", "headwind", "crosswind"]
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.7))
    x = np.arange(len(classes))
    width = 0.26
    for k, m in enumerate(METHODS[1:]):
        vals = [df[(df.wind_class == c) & (df.method == m)].extra_energy_pct.mean() for c in classes]
        axes[0].bar(x + (k - 1) * width, vals, width * 0.92, color=PALETTE[k + 1], label=METHOD_LABELS[m])
    axes[0].set_xticks(x, [c.capitalize() for c in classes])
    axes[0].set(ylabel="Extra energy vs ground truth (%)")
    axes[0].grid(axis="x", visible=False)
    axes[0].legend(fontsize=6.8)
    for k, m in enumerate(METHODS):
        e = np.sort(df[(df.method == m)].arrival_error_sec.abs().dropna().to_numpy())
        if len(e):
            axes[1].plot(e, np.linspace(0, 100, len(e)), color=PALETTE[k], label=METHOD_LABELS[m])
    axes[1].axvline(60, color="#52514e", ls="--", lw=0.8)
    axes[1].text(62, 4, "one 60-s slot", fontsize=7, color="#52514e")
    axes[1].set(xlabel="|Arrival - planned arrival| (s)", ylabel="Missions (%)")
    for k, m in enumerate(METHODS):
        vals = [df[(df.wind_class == c) & (df.method == m)].max_slot_shift.mean() for c in classes]
        axes[2].bar(x + (k - 1.5) * 0.2, vals, 0.18, color=PALETTE[k], label=METHOD_LABELS[m])
    axes[2].set_xticks(x, [c.capitalize() for c in classes])
    axes[2].set(ylabel="Reservation buffer needed (slots)")
    axes[2].grid(axis="x", visible=False)
    for ax, letter in zip(axes, "abc"):
        panel_label(ax, letter)
    return save(fig, out, f"value_of_wind_information_{tag}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--missions", type=int, default=256)
    ap.add_argument("--mission", type=int, default=None, help="run only mission k and print it")
    ap.add_argument("--wind-speed", type=float, default=15.0, help="reference wind at 400 ft (mph)")
    ap.add_argument("--select", default="energy", choices=["energy", "time", "balanced"],
                    help="rule picking one route from each Pareto set")
    ap.add_argument("--model", default="paper", choices=["paper", "continuous"],
                    help="wind projection: TRR piecewise rule or continuous along-track component")
    ap.add_argument("--lanes", type=int, default=2)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", default="results/wind_belief")
    a = ap.parse_args()

    corr = belief_corridor(n_lat=a.lanes)
    missions = sample_missions(corr, a.missions, a.seed)
    if a.mission is not None:
        m = missions[a.mission]
        rows, _, _, _ = run_mission(corr, VehicleModel(), m, a.wind_speed, a.select, a.model)
        print(f"Mission {m.index}: {m.wind_class}, wind toward {m.wind_dir:.0f} deg, heading {m.heading:.0f}")
        print(pd.DataFrame(rows)[["method", "planned_kwh", "realized_kwh", "extra_energy_pct",
                                  "planned_min", "realized_min", "max_slot_shift"]].round(4).to_string())
        return
    work = partial(_work, speed=a.wind_speed, rule=a.select, model=a.model, lanes=a.lanes)
    rows = []
    if a.workers > 1:
        with ProcessPoolExecutor(a.workers) as ex:
            for r in ex.map(work, missions, chunksize=4):
                rows += r
    else:
        for k, m in enumerate(missions):
            rows += work(m)
            if (k + 1) % 32 == 0:
                print(f"  {k + 1}/{len(missions)} missions", flush=True)
    df = pd.DataFrame(rows)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    tag = f"{a.wind_speed:g}mph_{a.select}"
    df.to_csv(out / f"missions_{tag}.csv", index=False)
    tab = summary_table(df)
    tab.to_csv(out / f"summary_{tag}.csv", index=False)
    figures(df, out / "figures", tag)
    pd.set_option("display.width", 160)
    print(tab.round(2).to_string(index=False))
    print(f"Results written to {out}/")


if __name__ == "__main__":
    main()

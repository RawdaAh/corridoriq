"""Paper 1 - wind-aware Pareto routing experiments (TRR, Sections 6.2-6.5).

    python run_wind_routing.py --experiment all
    python run_wind_routing.py --experiment orientation --phase-dependent

Experiments
  orientation : altitude-shear wind rotated in 45 deg steps          (Tables 1-2)
  spatial     : weak wind in blocks B1-B2, strong in B3-B5           (Tables 3-5)
  temporal    : calm->strong / strong->calm over a 30-min window     (Figure 7)
  congestion  : 0-90 % of directed links closed, 10 replications    (Figure 8)

Output: results/wind_routing/<experiment>.csv and figures/*.png|pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from aamsim.corridor import trr_wind_corridor
from aamsim.energy import VehicleModel, wind_link_costs
from aamsim.pareto import pareto_routes, route_profile
from aamsim.plotting import PALETTE, panel_label, plt, save
from aamsim.wind import trr_altitude_wind, trr_spatial_wind, trr_temporal_wind, zero_wind

DIRECTIONS = list(range(0, 360, 45))


def summarize(corr, routes):
    if not routes:
        return dict(n_paths=0, time_min=np.nan, energy_kwh=np.nan, ascents=np.nan, descents=np.nan)
    prof = [route_profile(corr, r) for r in routes]
    out = dict(n_paths=len(routes),
               time_min=np.mean([r.travel_min for r in routes]),
               energy_kwh=np.mean([r.energy_kwh for r in routes]),
               ascents=np.mean([p["ascents"] for p in prof]),
               descents=np.mean([p["descents"] for p in prof]))
    for key in ("alt_time_pct", "alt_energy_pct"):
        bands = {}
        for p in prof:
            for b, v in p[key].items():
                bands[b] = bands.get(b, 0) + v / len(prof)
        for b, v in bands.items():
            out[f"{key}[{b}]"] = v
    for key in ("block_time_pct", "block_energy_pct"):
        arr = np.mean([p[key] for p in prof], axis=0)
        for i, v in enumerate(arr, 1):
            out[f"{key}[B{i}]"] = v
    return out


def exp_orientation(corr, veh, o, d, spatial=False):
    rows = []
    maker = trr_spatial_wind if spatial else trr_altitude_wind
    cases = [("no_wind", zero_wind(corr))] + [(f"{a}", maker(corr, a)) for a in DIRECTIONS]
    for name, w in cases:
        rows.append(dict(wind_dir=name, **summarize(corr, pareto_routes(corr, veh, w, o, d))))
        print(f"  wind {name:>7}: {rows[-1]['n_paths']} paths, {rows[-1]['time_min']:.2f} min, "
              f"{rows[-1]['energy_kwh']:.4f} kWh")
    return pd.DataFrame(rows)


def exp_temporal(corr, veh, o, d, n_departures=10, bin_sec=180.0):
    rows = []
    for direction in (0, 180):
        for trend in ("calm_to_strong", "strong_to_calm"):
            w = trr_temporal_wind(corr, direction, trend, bin_sec=bin_sec)
            costs = wind_link_costs(corr, veh, w)
            for k in range(n_departures):
                R = pareto_routes(corr, veh, w, o, d, depart_sec=k * bin_sec, costs=costs)
                rows.append(dict(direction=direction, trend=trend, depart_min=k * bin_sec / 60,
                                 **summarize(corr, R)))
            print(f"  {direction:>3} deg {trend}: done")
    return pd.DataFrame(rows)


def exp_congestion(corr, veh, o, d, reps=10, seed=0):
    rows = []
    w = zero_wind(corr)
    costs = wind_link_costs(corr, veh, w)
    rng = np.random.default_rng(seed)
    for frac in np.round(np.arange(0, 1.0, 0.1), 1):
        for rep in range(reps):
            mask = rng.random(corr.n_arcs) >= frac
            R = pareto_routes(corr, veh, w, o, d, link_mask=mask, costs=costs)
            rows.append(dict(closed_frac=frac, rep=rep, feasible=len(R) > 0, **summarize(corr, R)))
        sub = [r for r in rows if r["closed_frac"] == frac]
        print(f"  closure {frac:.0%}: {sum(r['feasible'] for r in sub)}/{reps} feasible")
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------- figures

def fig_orientation(df, out, name):
    dfw = df[df.wind_dir != "no_wind"].copy()
    base = df[df.wind_dir == "no_wind"].iloc[0]
    labels = [{"0": "0\n(tail)", "180": "180\n(head)"}.get(x, x) for x in dfw.wind_dir]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6))
    for ax, col, lab, c, letter in [(axes[0], "time_min", "Average travel time (min)", PALETTE[0], "a"),
                                    (axes[1], "energy_kwh", "Average energy (kWh)", PALETTE[1], "b")]:
        ax.bar(range(len(dfw)), dfw[col], color=c, width=0.65)
        ax.axhline(base[col], color="#52514e", ls="--", lw=1)
        ax.text(len(dfw) - 0.5, base[col], " no wind", va="bottom", ha="right", fontsize=7.5,
                color="#52514e")
        ax.set_xticks(range(len(dfw)), labels, fontsize=7.5)
        ax.set(xlabel="Wind direction (deg)", ylabel=lab)
        ax.grid(axis="x", visible=False)
        panel_label(ax, letter)
    return save(fig, out, name)


def fig_temporal(df, out):
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8))
    styles = {(0, "calm_to_strong"): (PALETTE[0], "-"), (0, "strong_to_calm"): (PALETTE[2], "-"),
              (180, "calm_to_strong"): (PALETTE[7], "--"), (180, "strong_to_calm"): (PALETTE[1], "--")}
    for (dirn, trend), sub in df.groupby(["direction", "trend"]):
        c, ls = styles[(dirn, trend)]
        lab = f"{'Tail' if dirn == 0 else 'Head'}wind, {trend.replace('_', ' ')}"
        axes[0].plot(sub.depart_min, sub.time_min, color=c, ls=ls, marker="o", label=lab)
        axes[1].plot(sub.depart_min, sub.energy_kwh, color=c, ls=ls, marker="o", label=lab)
    axes[0].set(xlabel="Departure time (min)", ylabel="Average travel time (min)")
    axes[1].set(xlabel="Departure time (min)", ylabel="Average energy (kWh)")
    panel_label(axes[0], "a"); panel_label(axes[1], "b")
    axes[1].legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=7.5)
    return save(fig, out, "temporal_departure")


def fig_congestion(df, out):
    g = df.groupby("closed_frac")
    feas = g.feasible.mean()
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5))
    axes[0].bar(feas.index * 100, feas.values * 100, width=6, color=PALETTE[0])
    axes[0].set(xlabel="Links closed (%)", ylabel="Replications feasible (%)")
    for ax, col, lab, c in [(axes[1], "time_min", "Travel time (min)", PALETTE[0]),
                            (axes[2], "energy_kwh", "Energy (kWh)", PALETTE[1])]:
        m = g[col].mean()
        ax.plot(m.index * 100, m.values, marker="o", color=c)
        ax.set(xlabel="Links closed (%)", ylabel=lab)
    for ax, letter in zip(axes, "abc"):
        panel_label(ax, letter)
    return save(fig, out, "congestion")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--experiment", default="all",
                    choices=["orientation", "spatial", "temporal", "congestion", "all"])
    ap.add_argument("--phase-dependent", action="store_true",
                    help="add climb multipliers and the take-off energy term (TRR Sec. 5.4)")
    ap.add_argument("--lanes", type=int, default=2, help="lateral modules (lanes = value + 1)")
    ap.add_argument("--reps", type=int, default=10, help="replications for congestion")
    ap.add_argument("--out", default="results/wind_routing")
    a = ap.parse_args()

    corr = trr_wind_corridor(n_lat=a.lanes)
    veh = VehicleModel(phase_dependent=a.phase_dependent)
    o, d = corr.node(0, 0, 0), corr.node(corr.spec.n_long, 0, 0)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    tag = "_phase" if a.phase_dependent else ""
    print("Corridor:", corr.summary())
    todo = ["orientation", "spatial", "temporal", "congestion"] if a.experiment == "all" else [a.experiment]
    for e in todo:
        print(f"[{e}]")
        if e in ("orientation", "spatial"):
            df = exp_orientation(corr, veh, o, d, spatial=(e == "spatial"))
            df.to_csv(out / f"{e}{tag}.csv", index=False)
            fig_orientation(df, out / "figures", f"{e}{tag}")
        elif e == "temporal":
            df = exp_temporal(corr, veh, o, d)
            df.to_csv(out / f"temporal{tag}.csv", index=False)
            fig_temporal(df, out / "figures")
        else:
            df = exp_congestion(corr, veh, o, d, reps=a.reps)
            df.to_csv(out / f"congestion{tag}.csv", index=False)
            fig_congestion(df, out / "figures")
    print(f"Results written to {out}/")


if __name__ == "__main__":
    main()

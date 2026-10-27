"""Paper 2 - Shared-Intent Reservation Engine (SIRE) experiments (TRB paper).

    python run_sire_experiment.py --experiment demand --reps 30 --workers 8
    python run_sire_experiment.py --experiment all --reps 5        # quick look

Experiments (each varies one factor; common random numbers across cells)
  demand    : total demand 100 ... 4,000 trips                       (Figure 1)
  geometry  : lateral modules N_W 1-4 x altitude modules N_H 1-4     (Figure 2)
  direction : forward/reverse split x shared vs separated rule      (Figure 3)
  dispatch  : dispatch window 5-60 min x uniform vs peak pattern    (Figure 4)
  equity    : same-minute tie-break by operator id, lottery, rotation
  wind      : link times from a wind snapshot (tail vs head) - extension

Output: results/sire/<experiment>.csv (one row per replication) and figures/.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from aamsim.corridor import trb_sire_corridor
from aamsim.energy import VehicleModel, baseline_link_costs
from aamsim.plotting import PALETTE, band, panel_label, plt, save
from aamsim.sire import SireConfig, generate_trips, run_sire

STRESS = 2000
DEMANDS = [100, 200, 400, 600, 800, 1000, 1500, 2000, 3000, 4000]


def one_run(job):
    """Run one replication. job = dict(cell params). Returns a KPI row."""
    corr = trb_sire_corridor(job.get("n_long", 30), job.get("n_lat", 2), job.get("n_alt", 4))
    cfg = SireConfig(dispatch_window=job.get("window", 60),
                     conflict_rule=job.get("rule", "undirected"),
                     tie_break=job.get("tie_break", "operator"), seed=job["rep"])
    trips = generate_trips(corr, job["demand"], window=job.get("window", 60),
                           pattern=job.get("pattern", "uniform"),
                           forward_share=job.get("forward_share", 1.0),
                           seed=1000 + job["rep"])  # common random numbers across cells
    link_time = None
    if job.get("wind_dir") is not None:
        from aamsim.wind import trr_altitude_wind
        veh = VehicleModel()
        t0, _ = baseline_link_costs(corr, veh)
        lam, _ = trr_altitude_wind(corr, job["wind_dir"]).link_factors(corr, veh.cruise_mph)
        link_time = t0 * lam[0]
    res = run_sire(corr, trips, cfg, link_time_sec=link_time)
    return {**{k: v for k, v in job.items()}, **res.kpis()}


def run_jobs(jobs, workers):
    if workers <= 1:
        return [one_run(j) for j in jobs]
    with ProcessPoolExecutor(workers) as ex:
        return list(ex.map(one_run, jobs, chunksize=1))


def design(exp, reps):
    R = range(reps)
    if exp == "demand":
        return [dict(demand=d, rep=r) for d in DEMANDS for r in R]
    if exp == "geometry":
        return [dict(demand=STRESS, n_lat=w, n_alt=h, rep=r)
                for w in (1, 2, 3, 4) for h in (1, 2, 3, 4) for r in R]
    if exp == "direction":
        return [dict(demand=STRESS, forward_share=f, rule=rule, rep=r)
                for f in (1.0, 0.75, 0.5, 0.25) for rule in ("undirected", "directed") for r in R]
    if exp == "dispatch":
        return [dict(demand=STRESS, window=w, pattern=p, rep=r)
                for w in (5, 10, 20, 30, 45, 60) for p in ("uniform", "peak") for r in R]
    if exp == "equity":
        return [dict(demand=d, tie_break=tb, rep=r)
                for d in (1000, 2000) for tb in ("operator", "random", "rotating") for r in R]
    if exp == "wind":
        return [dict(demand=d, wind_dir=wd, rep=r)
                for d in (400, 1000, 2000) for wd in (None, 0, 180) for r in R]
    raise ValueError(exp)


# ------------------------------------------------------------------------- figures

def q05(x):
    return x.quantile(0.05)


def q95(x):
    return x.quantile(0.95)


def fig_demand(df, out):
    g = df.groupby("demand")
    x = g.success_rate.mean().index
    ref = df[df.demand == df.demand.min()].tt_acc_sec.mean()
    df = df.assign(itt=100 * (df.tt_acc_sec / ref - 1))
    g = df.groupby("demand")
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.0))
    for ax, col, lab, c, letter in [
            (axes[0, 0], "success_rate", "Success rate (SR)", PALETTE[0], "a"),
            (axes[0, 1], "cascade", "Cascade indicator (W = 5)", PALETTE[1], "b"),
            (axes[1, 0], "itt", "Travel-time inflation (%)", PALETTE[2], "c"),
            (axes[1, 1], "gini", "Gini of operator SR", PALETTE[6], "d")]:
        m = g[col].mean()
        band(ax, x, g[col].agg(q05), g[col].agg(q95), c)
        ax.plot(x, m, marker="o", color=c)
        ax.axvline(STRESS, color="#52514e", ls="--", lw=0.9)
        ax.set(xlabel="Total demand (trips)", ylabel=lab)
        panel_label(ax, letter)
    ax2 = axes[0, 0]
    acc = g.accepted.mean()
    for d in (STRESS, x.max()):
        if d in acc.index:
            ax2.annotate(f"{acc[d]:.0f} accepted", (d, g.success_rate.mean()[d]),
                         textcoords="offset points", xytext=(4, 9), fontsize=7, ha="left",
                         color="#52514e")
    return save(fig, out, "demand_response")


def fig_geometry(df, out):
    pv = df.groupby(["n_lat", "n_alt"]).success_rate.mean().unstack()
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.6))
    im = axes[0].imshow(pv.values * 100, origin="lower", cmap="Blues", aspect="auto",
                        extent=[0.5, 4.5, 0.5, 4.5], vmin=0, vmax=100)
    for i, w in enumerate(pv.index):
        for j, h in enumerate(pv.columns):
            v = pv.loc[w, h] * 100
            axes[0].text(h, w, f"{v:.0f}", ha="center", va="center", fontsize=7.5,
                         color="white" if v > 55 else "#0b0b0b")
    axes[0].set(xlabel="Altitude modules $N_H$", ylabel="Lateral modules $N_W$",
                xticks=[1, 2, 3, 4], yticks=[1, 2, 3, 4])
    axes[0].grid(False)
    fig.colorbar(im, ax=axes[0], fraction=0.046, label="SR (%)")
    for k, h in enumerate(pv.columns):
        d = np.diff(pv[h].values) * 100
        axes[1].plot(pv.index[1:], d, marker="o", color=PALETTE[k], label=f"$N_H$={h}")
    axes[1].set(xlabel="Lateral modules $N_W$ (to)", ylabel="SR gain per lateral module (pp)")
    for k, w in enumerate(pv.index):
        d = np.diff(pv.loc[w].values) * 100
        axes[2].plot(pv.columns[1:], d, marker="o", color=PALETTE[k], label=f"$N_W$={w}")
    axes[2].set(xlabel="Altitude modules $N_H$ (to)", ylabel="SR gain per altitude layer (pp)")
    lo = min(axes[1].get_ylim()[0], axes[2].get_ylim()[0])
    hi = max(axes[1].get_ylim()[1], axes[2].get_ylim()[1])
    for ax in axes[1:]:
        ax.set_ylim(lo, hi)
        ax.legend(fontsize=6.5)
    for ax, letter in zip(axes, "abc"):
        panel_label(ax, letter)
    return save(fig, out, "geometry")


def fig_direction(df, out):
    g = df.groupby(["forward_share", "rule"]).success_rate.mean().unstack() * 100
    g = g.sort_index(ascending=False)
    labels = [f"{int(f*100)}/{int(100-f*100)}" for f in g.index]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6))
    x = np.arange(len(g))
    axes[0].bar(x - 0.18, g["undirected"], 0.34, color=PALETTE[0], label="Shared (undirected)")
    axes[0].bar(x + 0.18, g["directed"], 0.34, color=PALETTE[1], label="Direction-separated")
    axes[0].set_xticks(x, labels)
    axes[0].set(xlabel="Forward / reverse split (%)", ylabel="Success rate (%)")
    axes[0].legend(fontsize=7)
    gain = g["directed"] - g["undirected"]
    axes[1].bar(x, gain, 0.5, color=PALETTE[2])
    for xi, v in zip(x, gain):
        axes[1].text(xi, v, f"{v:.1f}", ha="center", va="bottom", fontsize=7.5)
    axes[1].set_xticks(x, labels)
    axes[1].set(xlabel="Forward / reverse split (%)", ylabel="Gain from separation (pp)")
    for ax, letter in zip(axes, "ab"):
        ax.grid(axis="x", visible=False)
        panel_label(ax, letter)
    return save(fig, out, "direction")


def fig_dispatch(df, out):
    fig, ax = plt.subplots(figsize=(4.2, 2.8))
    for k, (p, sub) in enumerate(df.groupby("pattern")):
        g = sub.groupby("window").success_rate
        band(ax, g.mean().index, g.agg(q05) * 100, g.agg(q95) * 100, PALETTE[k])
        ax.plot(g.mean().index, g.mean() * 100, marker="o", color=PALETTE[k],
                ls="-" if p == "uniform" else "--", label=p.capitalize())
    ax.set(xlabel="Dispatch window (min)", ylabel="Success rate (%)")
    ax.legend()
    return save(fig, out, "dispatch_window")


def fig_equity(df, out):
    cols = [c for c in df.columns if c.startswith("sr_op")]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6))
    for ax, d, letter in zip(axes, sorted(df.demand.unique()), "ab"):
        sub = df[df.demand == d].groupby("tie_break")[cols].mean() * 100
        sub = sub.loc[[t for t in ("operator", "random", "rotating") if t in sub.index]]
        x = np.arange(len(sub))
        for k, c in enumerate(cols):
            ax.bar(x + (k - 1) * 0.25, sub[c], 0.23, color=PALETTE[k], label=f"Operator {c[-1]}")
        ax.set_xticks(x, ["By operator id", "Lottery", "Rotating"])
        ax.set(ylabel="Operator success rate (%)", title=f"Demand = {d:,} trips")
        ax.grid(axis="x", visible=False)
        panel_label(ax, letter)
    axes[1].legend(fontsize=7)
    return save(fig, out, "equity_tiebreak")


def fig_wind(df, out):
    df = df.assign(case=df.wind_dir.map({0: "Tailwind", 180: "Headwind"}).fillna("No wind"))
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6))
    for k, case in enumerate(["No wind", "Tailwind", "Headwind"]):
        g = df[df.case == case].groupby("demand")
        axes[0].plot(g.success_rate.mean().index, g.success_rate.mean() * 100, marker="o",
                     color=PALETTE[k], label=case)
        axes[1].plot(g.tt_acc_sec.mean().index, g.tt_acc_sec.mean() / 60, marker="o",
                     color=PALETTE[k], label=case)
    axes[0].set(xlabel="Total demand (trips)", ylabel="Success rate (%)")
    axes[1].set(xlabel="Total demand (trips)", ylabel="Accepted-trip travel time (min)")
    axes[1].legend(fontsize=7)
    for ax, letter in zip(axes, "ab"):
        panel_label(ax, letter)
    return save(fig, out, "wind_extension")


FIGS = dict(demand=fig_demand, geometry=fig_geometry, direction=fig_direction,
            dispatch=fig_dispatch, equity=fig_equity, wind=fig_wind)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--experiment", default="all", choices=list(FIGS) + ["all"])
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", default="results/sire")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    exps = list(FIGS) if a.experiment == "all" else [a.experiment]
    for e in exps:
        jobs = design(e, a.reps)
        print(f"[{e}] {len(jobs)} runs ...", flush=True)
        df = pd.DataFrame(run_jobs(jobs, a.workers))
        df.to_csv(out / f"{e}.csv", index=False)
        FIGS[e](df, out / "figures")
        keys = [k for k in ("demand", "n_lat", "n_alt", "forward_share", "rule", "window",
                            "pattern", "tie_break", "wind_dir") if k in df.columns and df[k].nunique() > 1]
        print(df.groupby(keys, dropna=False)[["success_rate", "accepted", "tt_acc_sec", "gini"]]
              .mean().round(3).to_string())
    print(f"Results written to {out}/")


if __name__ == "__main__":
    main()

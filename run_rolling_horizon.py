"""Paper 3 - priority-aware rolling-horizon corridor scheduling.

    python run_rolling_horizon.py --case base --seeds 3
    python run_rolling_horizon.py --design quick --seeds 5 --workers 4
    python run_rolling_horizon.py --design full --seeds 30 --workers 8   # notebook grid

Every case runs three methods on identical demand: causal FCFS, rolling
two-stage priority ILP, and the offline full-information restricted benchmark.
Every schedule is audited (geometry-rebuilt resource tokens, no conflicts,
request windows, freeze rule) before metrics are recorded.

Output: results/rolling_horizon/metrics.csv, summary.csv (paired ILP - FCFS
differences with bootstrap 95% CIs over seeds) and figures/.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from aamsim.plotting import PALETTE, panel_label, plt, save
from aamsim.rolling import RHConfig, experiment_design, run_case


def one(job):
    name, cfg, seed = job
    table, _, _, _ = run_case(replace(cfg, seed=seed))
    table["case"], table["seed"] = name, seed
    return table


def bootstrap_ci(x, seed=991, draws=5000):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        return (x.mean() if len(x) else np.nan, np.nan, np.nan)
    rng = np.random.default_rng(seed)
    m = rng.choice(x, size=(draws, len(x)), replace=True).mean(axis=1)
    return x.mean(), np.quantile(m, 0.025), np.quantile(m, 0.975)


def summarize(raw):
    rows = []
    metrics = ["weighted_service_rate", "service_rate", "mean_delay", "max_delay",
               "emergency_service_rate", "max_wall_sec"]
    for case, sub in raw.groupby("case", sort=False):
        piv = {m: sub.pivot(index="seed", columns="method", values=m) for m in metrics}
        for m in metrics:
            for a, b in (("ilp", "fcfs"), ("offline", "ilp")):
                if a in piv[m] and b in piv[m]:
                    mean, lo, hi = bootstrap_ci(piv[m][a] - piv[m][b])
                    rows.append(dict(case=case, metric=m, contrast=f"{a}_minus_{b}",
                                     mean=mean, ci_lo=lo, ci_hi=hi, n=len(piv[m])))
    return pd.DataFrame(rows)


def figures(raw, out):
    base = raw[raw.case.str.startswith("load_") | (raw.case == "base")].copy()
    base["n"] = base.requests
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.6))
    for k, m in enumerate(["fcfs", "ilp", "offline"]):
        g = base[base.method == m].groupby("n")
        lab = {"fcfs": "FCFS", "ilp": "Rolling ILP", "offline": "Offline benchmark"}[m]
        axes[0].plot(g.weighted_service_rate.mean().index, g.weighted_service_rate.mean() * 100,
                     marker="o", color=PALETTE[k], label=lab)
        axes[1].plot(g.mean_delay.mean().index, g.mean_delay.mean(), marker="o", color=PALETTE[k])
        axes[2].plot(g.emergency_service_rate.mean().index, g.emergency_service_rate.mean() * 100,
                     marker="o", color=PALETTE[k])
    axes[0].set(xlabel="Requests", ylabel="Priority-weighted service (%)")
    axes[1].set(xlabel="Requests", ylabel="Mean ground delay (min)")
    axes[2].set(xlabel="Requests", ylabel="Emergency service (%)")
    axes[0].legend(fontsize=7)
    for ax, letter in zip(axes, "abc"):
        panel_label(ax, letter)
    save(fig, out, "load_response")

    other = raw[~raw.case.str.startswith("load_")]
    piv = other.pivot_table(index="case", columns="method", values="weighted_service_rate")
    if {"ilp", "fcfs"} <= set(piv.columns) and len(piv) > 1:
        d = ((piv.ilp - piv.fcfs) * 100).sort_values()
        fig, ax = plt.subplots(figsize=(4.6, 0.25 * len(d) + 1.0))
        ax.barh(d.index, d.values, color=[PALETTE[0] if v >= 0 else PALETTE[1] for v in d.values])
        ax.axvline(0, color="#52514e", lw=0.8)
        ax.set(xlabel="ILP - FCFS weighted service (pp)")
        ax.grid(axis="y", visible=False)
        save(fig, out, "sensitivity")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--design", default="quick", choices=["base", "quick", "full"])
    ap.add_argument("--case", default=None, help="run a single named case from the design")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--first-seed", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", default="results/rolling_horizon")
    a = ap.parse_args()
    base = RHConfig()
    cases = [("base", base)] if a.design == "base" else experiment_design(base, quick=a.design == "quick")
    if a.case:
        cases = [c for c in cases if c[0] == a.case]
        if not cases:
            raise SystemExit(f"unknown case {a.case}")
    seeds = list(range(a.first_seed, a.first_seed + a.seeds))
    jobs = [(n, c, s) for n, c in cases for s in seeds]
    print(f"{len(cases)} cases x {len(seeds)} seeds x 3 methods = {3 * len(jobs)} method runs")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    tables = []
    if a.workers > 1:
        with ProcessPoolExecutor(a.workers) as ex:
            for t in ex.map(one, jobs):
                tables.append(t)
                print(f"  {t.case.iloc[0]} seed {t.seed.iloc[0]} done", flush=True)
    else:
        for j in jobs:
            tables.append(one(j))
            print(f"  {j[0]} seed {j[2]} done", flush=True)
    raw = pd.concat(tables, ignore_index=True)
    raw.to_csv(out / "metrics.csv", index=False)
    summ = summarize(raw)
    summ.to_csv(out / "summary.csv", index=False)
    figures(raw, out / "figures")
    print(raw.groupby(["case", "method"], sort=False)[
        ["weighted_service_rate", "service_rate", "mean_delay", "max_wall_sec"]].mean().round(3).to_string())
    print(f"Results written to {out}/")


if __name__ == "__main__":
    main()

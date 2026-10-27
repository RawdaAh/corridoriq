"""Print summary tables for any result CSV written by the run_*.py scripts.

    python summarize_results.py results/sire/demand.csv
    python summarize_results.py results/wind_belief/missions_15mph_energy.csv
    python summarize_results.py results/rolling_horizon/metrics.csv
"""
import sys

import pandas as pd

pd.set_option("display.width", 180)
pd.set_option("display.max_columns", 30)

KEYS = ["demand", "n_lat", "n_alt", "forward_share", "rule", "window", "pattern", "tie_break",
        "wind_dir", "case", "method", "wind_class", "closed_frac", "direction", "trend"]


def main(paths):
    for p in paths:
        df = pd.read_csv(p)
        print(f"\n=== {p}  ({len(df)} rows)")
        if "mission" in df.columns:          # wind-information experiment
            from run_wind_belief import summary_table
            print(summary_table(df).round(2).to_string(index=False))
            continue
        keys = [k for k in KEYS if k in df.columns and df[k].nunique(dropna=False) > 1]
        num = [c for c in df.select_dtypes("number").columns
               if c not in keys and c not in ("rep", "seed", "mission")]
        if keys:
            print(df.groupby(keys, dropna=False, sort=False)[num].mean().round(4).to_string())
        else:
            print(df[num].describe().round(4).to_string())


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])

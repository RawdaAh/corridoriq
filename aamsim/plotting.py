"""Shared figure style: paper-ready, colour-blind-validated categorical palette."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Fixed categorical order (validated for adjacent-pair CVD separation).
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"

plt.rcParams.update({
    "font.family": "serif", "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": INK2,
    "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.prop_cycle": matplotlib.cycler(color=PALETTE),
    "lines.linewidth": 1.8, "lines.markersize": 4.5, "legend.frameon": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300,
})


def save(fig, out_dir, name: str):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out / f"{name}.png", bbox_inches="tight")
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    return out / f"{name}.png"


def panel_label(ax, letter: str):
    ax.text(-0.14, 1.04, f"({letter})", transform=ax.transAxes, fontweight="bold", va="bottom")


def band(ax, x, lo, hi, color, alpha=0.18):
    ax.fill_between(x, lo, hi, color=color, alpha=alpha, linewidth=0)

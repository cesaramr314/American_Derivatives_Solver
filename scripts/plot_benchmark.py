"""Plots the output of amopt_benchmark.

Usage:
    python scripts/plot_benchmark.py results/benchmark.csv
Writes results/error_vs_time.png and results/iterations.png.
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

STYLE = {  # color + marker, so methods are not told apart by color alone
    "psor": {"color": "#2a78d6", "marker": "o", "label": "PSOR (ω = 1.5)"},
    "psor-adaptive": {"color": "#1f9e78", "marker": "^", "label": "PSOR (adaptive ω)"},
    "penalty": {"color": "#eb6834", "marker": "s", "label": "Penalty"},
}


def load(path):
    data = defaultdict(lambda: defaultdict(list))
    with open(path) as f:
        for row in csv.DictReader(f):
            m = row["method"]
            data[m]["n"].append(int(row["n_space"]))
            data[m]["err"].append(float(row["abs_error"]))
            data[m]["iters"].append(float(row["iters_per_step"]))
            data[m]["sec"].append(float(row["seconds"]))
    return data


def style_axes(ax):
    ax.grid(True, which="major", color="#e4e3dc", linewidth=0.8)
    ax.grid(False, which="minor")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(frameon=False)


def main():
    path = Path(sys.argv[1] if len(sys.argv) > 1 else "results/benchmark.csv")
    out = path.parent
    data = load(path)

    fig, ax = plt.subplots(figsize=(5.5, 4))
    for m, d in data.items():
        s = STYLE[m]
        ax.loglog(d["sec"], d["err"], color=s["color"], marker=s["marker"],
                  linewidth=2, markersize=7, label=s["label"])
    ax.set_xlabel("CPU time [s]")
    ax.set_ylabel("absolute error in V(S0)")
    ax.set_title("American put: error vs cost")
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(out / "error_vs_time.png", dpi=200)

    fig, ax = plt.subplots(figsize=(5.5, 4))
    for m, d in data.items():
        s = STYLE[m]
        ax.semilogx(d["n"], d["iters"], color=s["color"], marker=s["marker"],
                    linewidth=2, markersize=7, label=s["label"])
    ax.set_xlabel("space intervals N  (Δt ∝ Δx)")
    ax.set_ylabel("average iterations per time step")
    ax.set_title("LCP solver cost per time step")
    ax.set_ylim(bottom=0)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(out / "iterations.png", dpi=200)
    print(f"Figures saved in {out}/")


if __name__ == "__main__":
    main()
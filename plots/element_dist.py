"""Distribution of a chosen element's count per molecule, across datasets.

Ported from plot_functions.create_histogram_all + plot_kde_hist_all,
driven by ecomp.json files. Dropped the `mode`, `samples`, `valid` and
`compare` branches — current code paths only ever used mode='all'.

Two styles exposed to the UI:
  - 'histogram'      : grouped-ish density histograms (one per dataset)
  - 'histogram+kde'  : same histograms with a gaussian KDE overlaid
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from scipy.stats import gaussian_kde


def _element_counts(ecomp_path: str | Path, symbol: str) -> np.ndarray:
    """Count of `symbol` in every molecule (0 if absent).

    Mirrors plot_functions.get_element_counts: skips the trailing total
    entry, uses len(data)-1 as the molecule count.
    """
    data = json.loads(Path(ecomp_path).read_text())
    n = len(data) - 1
    out = np.zeros(n, dtype=int)
    for i in range(n):
        out[i] = data[i].get(symbol, 0)
    return out


def _available_elements(series: list[tuple[str, str]]) -> list[str]:
    """Union of every element that appears in any of the chosen datasets."""
    syms: set[str] = set()
    for _, path in series:
        data = json.loads(Path(path).read_text())
        for mol in data[:-1]:
            syms.update(mol)
    return sorted(syms)


def element_distribution(
    series: list[tuple[str, str]],
    symbol: str,
    style: str = "histogram",
) -> Figure:
    """series: [(label, ecomp_json_path), ...]
    symbol: element to plot
    style:  'histogram' or 'histogram+kde'
    """
    if not series:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.text(0.5, 0.5, "Pick at least one dataset", ha="center", va="center")
        ax.set_axis_off()
        return fig

    arrays: list[np.ndarray] = [_element_counts(p, symbol) for _, p in series]
    max_count = max((a.max() if a.size else 0) for a in arrays)

    fig, ax = plt.subplots(figsize=(max(8.0, 0.5 * (max_count + 2) + 2.0), 5.0))

    # No molecule contains this element in any picked dataset — render a
    # placeholder instead of letting matplotlib collapse to zero height.
    if max_count == 0:
        ax.text(
            0.5, 0.5,
            f"No '{symbol}' atoms in any picked dataset.",
            ha="center", va="center", fontsize=14,
            transform=ax.transAxes,
        )
        ax.set_axis_off()
        fig.tight_layout()
        return fig

    # Integer-centred bins so each count sits in its own bar.
    bins = np.arange(max_count + 2) - 0.5
    ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.45, zorder=0)
    ax.set_axisbelow(True)

    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", ["#4C72B0"])
    for i, ((label, _), arr) in enumerate(zip(series, arrays)):
        color = cycle[i % len(cycle)]
        ax.hist(
            arr, bins=bins, density=True,
            histtype="stepfilled",
            edgecolor="black", linewidth=0.6,
            color=color, alpha=0.45,
            label=label, zorder=3,
        )
        if style == "histogram+kde" and arr.size > 1 and arr.max() > 0:
            try:
                density = gaussian_kde(arr, bw_method=0.3)
                x = np.linspace(0, max_count, 300)
                ax.plot(x, density(x), color=color, linewidth=2, zorder=4)
            except np.linalg.LinAlgError:
                # KDE fails if the data is constant — just skip the overlay.
                pass

    ax.set_xlabel(f"{symbol} atoms per molecule", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_xticks(range(0, max_count + 1))
    ax.tick_params(axis="both", labelsize=10)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    if len(series) > 4:
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
                  frameon=False, fontsize=10)
        fig.tight_layout(rect=[0, 0, 0.86, 1])
    else:
        ax.legend(frameon=False, fontsize=10)
        fig.tight_layout()
    return fig

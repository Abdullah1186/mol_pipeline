"""Average-atom-proportion bar plot, driven by ecomp.json files.

A port of `create_bar_all` from plot_functions.py, cleaned up for
Streamlit use:
  - takes explicit (label, json_path) pairs instead of a hardcoded
    base_dir + tuple-of-tuples
  - returns a matplotlib Figure (no plt.show())
  - dropped 'compare' / 'split' / 'samples' modes — only 'mode=all' is
    used and the dead branches were confusing

Data shape expected in each ecomp.json:
  list of {symbol: count} dicts, with the LAST entry being the running
  total over the whole set. This matches what pipelines/steps/write_jsons.py
  produces and what plot_functions.generate_data wrote historically.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ase.data import atomic_numbers, chemical_symbols
from matplotlib.figure import Figure


def _average_composition(ecomp_path: str | Path) -> dict[str, float]:
    """Per-molecule average elemental proportion. Mirrors
    plot_functions.get_average_composition: per molecule we compute the
    fraction (count / total atoms in that molecule) for each element,
    sum across molecules, then divide by molecule count.

    The last entry in the JSON list is the dataset total — we skip it for
    the per-molecule averaging but use len(list) - 1 as the denominator
    (matches the original)."""
    data = json.loads(Path(ecomp_path).read_text())
    avg: dict[str, float] = {}
    for molecule in data[:-1]:  # exclude the trailing total
        n_atoms = sum(molecule.values())
        if n_atoms == 0:
            continue
        for sym, count in molecule.items():
            avg[sym] = avg.get(sym, 0.0) + count / n_atoms
    denom = len(data) - 1
    if denom > 0:
        for sym in avg:
            avg[sym] /= denom
    return avg


def _average_composition_by_atomic_number(ecomp_path: str | Path) -> dict[int, float]:
    """Same as above but rekeyed by atomic number (so the bar plot can
    sort along the x-axis by Z)."""
    by_sym = _average_composition(ecomp_path)
    return {atomic_numbers[sym]: frac for sym, frac in by_sym.items()}


def bar_average_proportion(series: list[tuple[str, str]]) -> Figure:
    """Grouped bar plot of average atom proportion per molecule.

    series: list of (label, ecomp_json_path). One group of bars per series.
    Returns: matplotlib Figure.

    Math identical to the original create_bar_all (via _average_composition);
    visuals tuned for Streamlit: clean grid, percent y-axis, symbol+Z xticks,
    optional value labels when the plot is sparse enough to fit them.
    """
    if not series:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.text(0.5, 0.5, "Pick at least one dataset", ha="center", va="center")
        ax.set_axis_off()
        return fig

    # Compute the per-series {Z: avg_frac} and pool all Zs that appear.
    avgs: list[dict[int, float]] = [
        _average_composition_by_atomic_number(p) for _, p in series
    ]
    all_z = sorted({z for a in avgs for z in a})

    # Pad zeros where an element is absent from a series.
    for a in avgs:
        for z in all_z:
            a.setdefault(z, 0.0)

    # Matrix shape (n_series, n_elements), sorted by Z. Multiply by 100 so
    # the y-axis reads as % (math is the same — purely a display rescale).
    pct = np.array([[a[z] * 100 for z in all_z] for a in avgs])

    n_series = len(series)
    n_el = len(all_z)
    x = np.arange(n_el)
    bar_width = min(0.85 / n_series, 0.28)
    shifts = (np.arange(n_series) - (n_series - 1) / 2) * bar_width

    # Width scales with element count; height stays fixed.
    fig_w = max(7.0, 0.9 * n_el + 1.5 + 0.2 * n_series)
    fig, ax = plt.subplots(figsize=(fig_w, 4.8))

    # Light horizontal grid only — clean look.
    ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.45, zorder=0)
    ax.set_axisbelow(True)

    # Use the default mpl color cycle but with a touch of transparency and a
    # thin edge so adjacent bars stay distinguishable.
    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", ["#4C72B0"])
    show_values = n_series <= 3 and n_el <= 12

    for i, (label, _) in enumerate(series):
        color = cycle[i % len(cycle)]
        bars = ax.bar(
            x + shifts[i], pct[i], bar_width,
            label=label,
            color=color,
            edgecolor="black",
            linewidth=0.6,
            alpha=0.88,
            zorder=3,
        )
        if show_values:
            for bar, v in zip(bars, pct[i]):
                if v < 0.5:  # skip near-zero labels to avoid clutter
                    continue
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    v,
                    f"{v:.1f}",
                    ha="center", va="bottom",
                    fontsize=8, color="#333",
                )

    ax.set_xlabel("Element", fontsize=12)
    ax.set_ylabel("Average proportion per molecule (%)", fontsize=12)
    ax.set_xticks(x)
    # Two-line tick label: element symbol on top, atomic number below.
    labels = [f"{chemical_symbols[z]}\n$Z={z}$" for z in all_z]
    ax.set_xticklabels(labels, fontsize=11)
    ax.tick_params(axis="y", labelsize=10)
    ax.set_ylim(0, max(pct.max() * 1.18, 1))
    ax.margins(x=0.02)

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    # Legend outside on the right when there are many series; inline otherwise.
    if n_series > 4:
        ax.legend(
            loc="center left", bbox_to_anchor=(1.01, 0.5),
            frameon=False, fontsize=10,
        )
        fig.tight_layout(rect=[0, 0, 0.86, 1])
    else:
        ax.legend(frameon=False, fontsize=10, loc="upper right")
        fig.tight_layout()
    return fig

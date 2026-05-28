"""Plot registry.

Every plot is a function:

    render(named_datasets: dict[str, list[(positions, atom_types)]],
           dataset_info: dict) -> matplotlib.figure.Figure | list[Figure]

The dict keys are user-facing labels (the `name` field from params.yaml /
the Streamlit UI) so the figure can show one series per DB. Single
matplotlib axes overlaid by default; functions are free to return a list
of figures if overlay isn't meaningful for them.

To add a real plot:
  1. Write a function with the signature above (or import an existing one
     and write a thin adapter).
  2. Add it to PLOTS below with a human-readable key.
"""

from __future__ import annotations

from collections import Counter
from typing import Callable

import matplotlib.pyplot as plt
from matplotlib.figure import Figure


def _atom_count_per_molecule(
    named_datasets: dict[str, list], dataset_info: dict
) -> Figure:
    """Placeholder plot: histogram of #atoms per molecule, overlaid per DB.

    Cheap to compute, exercises the overlay contract end-to-end so the UI
    is testable before real plots are wired in.
    """
    fig, ax = plt.subplots(figsize=(7, 4))
    for name, mols in named_datasets.items():
        counts = [int(at.shape[0]) for _, at in mols]
        ax.hist(counts, bins=30, alpha=0.5, label=f"{name}  (n={len(mols)})")
    ax.set_xlabel("atoms per molecule")
    ax.set_ylabel("count")
    ax.set_title("Atom count distribution")
    ax.legend()
    fig.tight_layout()
    return fig


def _element_distribution(
    named_datasets: dict[str, list], dataset_info: dict
) -> Figure:
    """Placeholder plot: total atom counts per element symbol, grouped bars."""
    decoder = dataset_info["atom_decoder"]
    fig, ax = plt.subplots(figsize=(7, 4))

    # Tally per-element counts for each dataset.
    series = {}
    all_symbols: set[str] = set()
    for name, mols in named_datasets.items():
        c: Counter = Counter()
        for _, at in mols:
            for a in at.tolist():
                c[decoder[int(a)]] += 1
        series[name] = c
        all_symbols |= set(c)

    symbols = sorted(all_symbols)
    width = 0.8 / max(len(series), 1)
    for i, (name, c) in enumerate(series.items()):
        ax.bar(
            [j + i * width for j in range(len(symbols))],
            [c.get(s, 0) for s in symbols],
            width=width,
            label=name,
        )
    ax.set_xticks([j + 0.4 - width / 2 for j in range(len(symbols))])
    ax.set_xticklabels(symbols)
    ax.set_xlabel("element")
    ax.set_ylabel("atom count")
    ax.set_title("Element distribution")
    ax.legend()
    fig.tight_layout()
    return fig


# The registry the UI iterates. Edit this one line per new plot.
PLOTS: dict[str, Callable] = {
    "Atom count per molecule": _atom_count_per_molecule,
    "Element distribution": _element_distribution,
}

"""Distributions of atom count and molecular weight per molecule.

Ported from plot_functions.plot_atoms_data_all + get_atoms_data.
Returns two figures (one per quantity). All data comes from ecomp.json;
no DB access needed. Dropped `valid`, `samples`, `split`, `homo_lumo` —
current paths only use the plain 'all' version.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ase.data import atomic_masses, atomic_numbers
from matplotlib.figure import Figure
from scipy.stats import gaussian_kde


def _atom_counts_and_weights(ecomp_path: str | Path) -> tuple[list[int], list[float]]:
    """Per-molecule atom count and molecular weight. Mirrors
    get_atoms_data without the valid/sample/split/homo_lumo branches."""
    data = json.loads(Path(ecomp_path).read_text())
    no_of_atoms = [sum(mol.values()) for mol in data[:-1]]
    weights = [
        sum(atomic_masses[atomic_numbers[sym]] * n for sym, n in mol.items())
        for mol in data[:-1]
    ]
    return no_of_atoms, weights


def _kde_figure(
    series_values: list[tuple[str, list[float]]],
    xlabel: str,
    xlim: tuple[float, float] | None = None,
) -> Figure:
    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.45, zorder=0)
    ax.set_axisbelow(True)

    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", ["#4C72B0"])
    drew_any = False
    for i, (label, vals) in enumerate(series_values):
        if len(vals) < 2 or max(vals) == min(vals):
            continue  # KDE requires variance
        color = cycle[i % len(cycle)]
        try:
            density = gaussian_kde(vals, bw_method=0.3)
        except np.linalg.LinAlgError:
            continue
        x = np.linspace(0, max(vals), 1000)
        ax.plot(x, density(x), label=label, color=color, alpha=0.85, linewidth=2)
        drew_any = True

    if not drew_any:
        # Every series was skipped (constant data, too few points, or empty).
        # Draw a message so the tile keeps a reasonable size.
        ax.clear()
        ax.text(
            0.5, 0.5,
            "No plottable data (need ≥2 points with variance per dataset).",
            ha="center", va="center", fontsize=13,
            transform=ax.transAxes,
        )
        ax.set_axis_off()
        fig.tight_layout()
        return fig

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.tick_params(axis="both", labelsize=10)
    if xlim:
        ax.set_xlim(*xlim)
    else:
        ax.margins(x=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    if len(series_values) > 4:
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
                  frameon=False, fontsize=10)
        fig.tight_layout(rect=[0, 0, 0.86, 1])
    else:
        ax.legend(frameon=False, fontsize=10)
        fig.tight_layout()
    return fig


def atom_count_distribution(series: list[tuple[str, str]]) -> Figure:
    """KDE of #atoms per molecule, one curve per dataset."""
    data = [(label, _atom_counts_and_weights(p)[0]) for label, p in series]
    return _kde_figure(data, xlabel="Atoms per molecule")


def weight_distribution(series: list[tuple[str, str]]) -> Figure:
    """KDE of molecular weight (g/mol), one curve per dataset."""
    data = [(label, _atom_counts_and_weights(p)[1]) for label, p in series]
    return _kde_figure(data, xlabel="Molecular weight (g/mol)")

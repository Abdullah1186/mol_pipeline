"""Evaluate: V/U/N measurement + per-molecule filter flags.

Splits the old keep-set logic in two:
  - measurement (validity / uniqueness / novelty) stays here and always runs
  - per-molecule masks (`is_valid`, `is_unique`) are emitted as outputs so
    a downstream `ApplyFilters` step can AND user-selected criteria

This step never drops anything itself. That's `ApplyFilters`' job.
"""

from __future__ import annotations

from typing import Any

from rdkit import Chem

from rdkit_functions import BasicMolecularMetrics, build_molecule, mol2smiles

from ..context import RunContext
from ..datasets import info_for
from .base import Step


class Evaluate(Step):
    name = "evaluate"
    inputs = ("molecules", "source_ids")
    outputs = (
        "validity", "uniqueness", "novelty",
        "valid_count", "unique_count", "novel_count",
        "is_valid", "is_unique",
    )

    def run(self, ctx: RunContext, *, molecules, source_ids, **_: Any) -> dict[str, Any]:
        info = info_for(ctx.input.dataset)

        # Largest-fragment SMILES per molecule (None if invalid).
        db_smiles: list[str | None] = []
        for positions, atom_types in molecules:
            mol = build_molecule(positions, atom_types, info)
            smi = mol2smiles(mol)
            if smi is not None:
                frags = Chem.rdmolops.GetMolFrags(mol, asMols=True)
                largest = max(frags, default=mol, key=lambda m: m.GetNumAtoms())
                smi = mol2smiles(largest)
            db_smiles.append(smi)

        metrics = BasicMolecularMetrics(info)
        (validity, uniqueness, novelty), _, _, stats = metrics.evaluate(
            molecules, check_odd_e=False
        )

        is_valid: list[bool] = [s is not None for s in db_smiles]

        # Uniqueness mask: mark the LAST index per SMILES as the survivor
        # (matches the dict-comprehension semantics in old filter.py:109-113,
        # which the existing baselines were captured against).
        last_index_per_smiles: dict[str, int] = {}
        for i, s in enumerate(db_smiles):
            if s is not None:
                last_index_per_smiles[s] = i
        survivors = set(last_index_per_smiles.values())
        is_unique: list[bool] = [i in survivors for i in range(len(molecules))]

        return {
            "validity": validity,
            "uniqueness": uniqueness,
            "novelty": novelty,
            "valid_count": stats.get("valid", 0),
            "unique_count": stats.get("unique", 0),
            "novel_count": stats.get("novel", 0),
            "is_valid": is_valid,
            "is_unique": is_unique,
        }

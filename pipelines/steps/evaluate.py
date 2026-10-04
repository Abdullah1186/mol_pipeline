"""Evaluate: V/U/N measurement + per-molecule filter flags.

Splits the old keep-set logic in two:
  - measurement (validity / uniqueness / novelty) stays here and always runs
  - per-molecule masks (`is_valid`, `is_unique`) are emitted as outputs so
    a downstream `ApplyFilters` step can AND user-selected criteria

Also emits `db_smiles` — the largest-fragment SMILES for every molecule
(or None if invalid) — so downstream WriteJSONs steps can reuse them
instead of recomputing (big win, RDKit was the dominant pipeline cost).

This step never drops anything itself. That's `ApplyFilters`' job.
"""

from __future__ import annotations

from typing import Any

from rdkit import Chem

from rdkit_functions import build_molecule, mol2smiles

from ..context import RunContext
from ..datasets import info_for
from .base import Step


class Evaluate(Step):
    name = "evaluate"
    inputs = ("molecules", "source_ids")
    outputs = (
        "validity", "uniqueness", "novelty",
        "valid_count", "unique_count", "novel_count",
        "is_valid", "is_unique", "db_smiles",
    )

    def run(self, ctx: RunContext, *, molecules, source_ids, **_: Any) -> dict[str, Any]:
        info = info_for(ctx.input.dataset)
        n = len(molecules)

        # Single RDKit pass per molecule: build mol, sanitize to SMILES,
        # then re-smiles over the largest fragment. db_smiles[i] is None iff
        # the molecule is invalid. This replaces TWO independent RDKit passes
        # the old code made (one here, one inside BasicMolecularMetrics.evaluate).
        db_smiles: list[str | None] = []
        for positions, atom_types in molecules:
            mol = build_molecule(positions, atom_types, info)
            smi = mol2smiles(mol)
            if smi is not None:
                frags = Chem.rdmolops.GetMolFrags(mol, asMols=True)
                largest = max(frags, default=mol, key=lambda m: m.GetNumAtoms())
                smi = mol2smiles(largest)
            db_smiles.append(smi)

        is_valid: list[bool] = [s is not None for s in db_smiles]
        valid_count = sum(is_valid)
        validity = valid_count / n if n else 0.0

        # Uniqueness mask: mark the LAST index per SMILES as the survivor
        # (matches the dict-comprehension semantics in old filter.py:109-113,
        # which the baselines were captured against).
        last_index_per_smiles: dict[str, int] = {}
        for i, s in enumerate(db_smiles):
            if s is not None:
                last_index_per_smiles[s] = i
        survivors = set(last_index_per_smiles.values())
        is_unique: list[bool] = [i in survivors for i in range(n)]
        unique_count = len(survivors)
        uniqueness = unique_count / valid_count if valid_count else 0.0

        # Novelty is not meaningful here (no reference SMILES set loaded) —
        # keep the fields so write_metrics doesn't need to change.
        novel_count = 0
        novelty = 0.0

        # Match the print lines the old BasicMolecularMetrics.evaluate emitted
        # so terminal output stays recognisable.
        print(f"Validity over {n} molecules: {validity * 100:.2f}%")
        if valid_count > 0:
            print(f"duplicates: {valid_count - unique_count}")
            print(f"Uniqueness over {valid_count} valid molecules: {uniqueness * 100:.2f}%")

        return {
            "validity": validity,
            "uniqueness": uniqueness,
            "novelty": novelty,
            "valid_count": valid_count,
            "unique_count": unique_count,
            "novel_count": novel_count,
            "is_valid": is_valid,
            "is_unique": is_unique,
            "db_smiles": db_smiles,
        }

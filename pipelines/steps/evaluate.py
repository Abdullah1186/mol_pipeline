"""Evaluate: validity, uniqueness, novelty + which molecules survive.

Wraps the body of filter.py:87-117. Does two things the old code did
inline:
  1. Calls BasicMolecularMetrics.evaluate() for the V/U/N metrics.
  2. Reproduces filter.py's DB_smiles -> keep_index logic so downstream
     steps get a smaller `kept_molecules` / `kept_source_ids` list.

We keep the original DB_smiles build (largest-fragment SMILES) even
though BasicMolecularMetrics.evaluate computes SMILES internally — the
old code did both, and the keep_index mapping uses the largest-frag
version. Behavior must match.
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
        "kept_molecules", "kept_source_ids",
    )

    def run(self, ctx: RunContext, *, molecules, source_ids, **_: Any) -> dict[str, Any]:
        info = info_for(ctx.input.dataset)

        # Largest-fragment SMILES per molecule (filter.py:88-96)
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
        (validity, uniqueness, novelty), unique_smiles, _novel_smiles, stats = metrics.evaluate(
            molecules, check_odd_e=False
        )

        # Map unique SMILES back to molecule indices (filter.py:109-117).
        # Note: dict-comprehension dedup keeps the LAST index per smile,
        # matching filter.py exactly.
        kept_mols = []
        kept_ids = []
        if unique_smiles:
            smiles_to_index = {s: i for i, s in enumerate(db_smiles) if s is not None}
            keep_index = {smiles_to_index[s] for s in unique_smiles if s in smiles_to_index}
            for idx in keep_index:
                kept_mols.append(molecules[idx])
                kept_ids.append(source_ids[idx])

        return {
            "validity": validity,
            "uniqueness": uniqueness,
            "novelty": novelty,
            "valid_count": stats.get("valid", 0),
            "unique_count": stats.get("unique", 0),
            "novel_count": stats.get("novel", 0),
            "kept_molecules": kept_mols,
            "kept_source_ids": kept_ids,
        }

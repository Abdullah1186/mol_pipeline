"""WriteJSONs: emit smiles.json and ecomp.json for a molecule set.

Two artifacts per call:
  - <basename>.smiles.json : list of SMILES strings (one per molecule).
                              SMILES is None when RDKit can't sanitise the
                              molecule — stored as JSON null.
  - <basename>.ecomp.json  : list of {element: count} dicts (one per
                              molecule) PLUS a final dict that is the
                              total over the whole set. Mirrors the layout
                              of plot_functions.generate_data so existing
                              plotting code that reads these works unchanged.

The runner calls this twice per input — once for the raw load and once
for the post-filter kept set — with different basenames so they don't
collide.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from rdkit_functions import build_molecule, mol2smiles

from ..context import RunContext
from ..datasets import info_for
from .base import Step


class WriteJSONs(Step):
    inputs = ("molecules",)
    outputs = ("smiles_path", "ecomp_path")

    def __init__(self, tag: str, basename: str):
        # tag: 'raw' or 'filtered' — only affects the step name in logs/manifest
        # basename: the on-disk stem; outputs become <basename>.smiles.json etc.
        self.name = f"write_jsons_{tag}"
        self._basename = basename

    def run(self, ctx: RunContext, *, molecules, **_: Any) -> dict[str, Any]:
        info = info_for(ctx.input.dataset)
        decoder = info["atom_decoder"]
        out_dir = Path(ctx.params.output_dir)

        # SMILES per molecule (None if invalid). Matches plot_functions.mol2smiles
        # semantics, including the 'largest fragment' rule used elsewhere in the
        # pipeline (see evaluate.py) — but here we want the unmodified SMILES
        # for the raw molecule, since that's what generate_data wrote. Keeping
        # it simple: one mol2smiles call, no largest-fragment selection.
        smiles_list: list[str | None] = []
        for positions, atom_types in molecules:
            mol = build_molecule(positions, atom_types, info)
            smiles_list.append(mol2smiles(mol))

        # Per-molecule element counts + a running total.
        all_ecomp: list[dict[str, int]] = []
        total: dict[str, int] = {}
        for _, atom_types in molecules:
            counts: dict[str, int] = {}
            for a in atom_types.tolist():
                sym = decoder[int(a)]
                counts[sym] = counts.get(sym, 0) + 1
                total[sym] = total.get(sym, 0) + 1
            all_ecomp.append(counts)
        all_ecomp.append(total)

        smiles_path = out_dir / f"{self._basename}.smiles.json"
        ecomp_path = out_dir / f"{self._basename}.ecomp.json"
        smiles_path.write_text(json.dumps(smiles_list))
        ecomp_path.write_text(json.dumps(all_ecomp))
        return {"smiles_path": str(smiles_path), "ecomp_path": str(ecomp_path)}

"""LoadDB: read the source ASE .db into in-memory molecules."""

from __future__ import annotations

from typing import Any

import torch
from ase.db import connect

from ..context import RunContext
from ..datasets import atom_encoder_for
from .base import Step


class LoadDB(Step):
    name = "load_db"
    inputs = ()
    outputs = ("molecules", "source_ids")

    def run(self, ctx: RunContext, **_: Any) -> dict[str, Any]:
        encoder = atom_encoder_for(ctx.input.dataset)

        molecules = []
        source_ids = []
        for row in connect(ctx.input.path).select():
            atoms = row.toatoms()
            positions = torch.tensor(atoms.positions, dtype=torch.float32)
            atom_types = torch.tensor([encoder[s] for s in atoms.get_chemical_symbols()])
            molecules.append((positions, atom_types))
            source_ids.append(row.id)

        return {"molecules": molecules, "source_ids": source_ids}

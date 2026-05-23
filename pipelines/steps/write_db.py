"""WriteDB: write kept molecules to a new ASE .db file.

Reads the source DB by path to fetch original Atoms objects (preserves
any metadata our (positions, atom_types) tuple discarded).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ase.db import connect

from .. import catalog
from ..context import RunContext
from .base import Step


class WriteDB(Step):
    name = "write_db"
    inputs = ("kept_source_ids",)
    outputs = ("output_path",)

    def run(self, ctx: RunContext, *, kept_source_ids, **_: Any) -> dict[str, Any]:
        if not kept_source_ids:
            return {"output_path": None}

        out_name = catalog.filtered_db_name(
            ctx.input.name, filter_odd_e=ctx.params.filter_odd_e
        )
        out_path = Path(ctx.params.output_dir) / out_name

        src = connect(ctx.input.path)
        dst = connect(str(out_path))
        for sid in kept_source_ids:
            dst.write(src.get(sid).toatoms())

        return {"output_path": str(out_path)}

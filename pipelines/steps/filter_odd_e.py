"""FilterOddE: drop odd-electron molecules from the kept set.

No-op when params.filter_odd_e is false (returns inputs unchanged).
Wraps filter.py:125-131.
"""

from __future__ import annotations

from typing import Any

from ..context import RunContext
from ..datasets import has_odd_electrons, info_for
from .base import Step


class FilterOddE(Step):
    name = "filter_odd_e"
    inputs = ("kept_molecules", "kept_source_ids")
    outputs = ("kept_molecules", "kept_source_ids")

    def run(
        self, ctx: RunContext, *, kept_molecules, kept_source_ids, **_: Any
    ) -> dict[str, Any]:
        if not ctx.params.filter_odd_e:
            return {"kept_molecules": kept_molecules, "kept_source_ids": kept_source_ids}

        info = info_for(ctx.input.dataset)
        mask = [not has_odd_electrons(at, info) for _, at in kept_molecules]
        return {
            "kept_molecules": [m for m, k in zip(kept_molecules, mask) if k],
            "kept_source_ids": [s for s, k in zip(kept_source_ids, mask) if k],
        }

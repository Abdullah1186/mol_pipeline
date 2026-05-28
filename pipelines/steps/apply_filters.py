"""ApplyFilters: AND together user-selected criteria into one mask.

Inputs: molecules + source_ids + per-mol bool masks from Evaluate.
Output: kept_molecules + kept_source_ids.

Each criterion is governed by a flag on ctx.params.filters. Criteria
not selected are treated as "all True" (i.e. don't constrain).
"""

from __future__ import annotations

from typing import Any

from ..context import RunContext
from ..datasets import has_odd_electrons, info_for
from .base import Step


class ApplyFilters(Step):
    name = "apply_filters"
    inputs = ("molecules", "source_ids", "is_valid", "is_unique")
    outputs = ("kept_molecules", "kept_source_ids")

    def run(
        self, ctx: RunContext, *, molecules, source_ids, is_valid, is_unique, **_: Any
    ) -> dict[str, Any]:
        flags = ctx.params.filters
        info = info_for(ctx.input.dataset)
        n = len(molecules)

        keep_valid = is_valid if flags.valid else [True] * n
        keep_unique = is_unique if flags.unique else [True] * n
        keep_even = (
            [not has_odd_electrons(at, info) for _, at in molecules]
            if flags.even_electrons else [True] * n
        )

        mask = [v and u and e for v, u, e in zip(keep_valid, keep_unique, keep_even)]
        kept_mols = [m for m, k in zip(molecules, mask) if k]
        kept_ids = [s for s, k in zip(source_ids, mask) if k]
        return {"kept_molecules": kept_mols, "kept_source_ids": kept_ids}

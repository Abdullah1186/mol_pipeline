"""ScanOddE: count molecules whose total electron count is odd.

Pure stat — does not mutate the molecule list. Used twice in the
pipeline: once on the raw load (Odd_e_before_filtering_VUN in the CSV)
and once on the V/U/N-kept set (Odd_e_after_filtering_VUN). The runner
distinguishes the two by passing a `tag` kwarg into the step instance.
"""

from __future__ import annotations

from typing import Any

from ..context import RunContext
from ..datasets import has_odd_electrons, info_for
from .base import Step


class ScanOddE(Step):
    inputs = ("molecules",)
    outputs = ("odd_count",)

    def __init__(self, tag: str):
        # tag = 'before' or 'after' — only affects the step name in logs/manifest
        self.name = f"scan_odd_e_{tag}"

    def run(self, ctx: RunContext, *, molecules, **_: Any) -> dict[str, Any]:
        info = info_for(ctx.input.dataset)
        count = sum(1 for _, at in molecules if has_odd_electrons(at, info))
        return {"odd_count": count}

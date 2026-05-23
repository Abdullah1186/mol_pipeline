"""Step ABC: the contract every pipeline step honors.

A Step:
  - has a name (used in logs and the manifest)
  - declares its input and output artifact names
  - implements run(ctx, **inputs) -> dict[str, Any]

The runner uses `inputs`/`outputs` to (a) verify wiring before executing
anything (catch typos early) and (b) thread artifacts between steps by
name. This is the same idea Airflow/Prefect/Dagster encode with task
decorators — we're just doing it by hand in 30 lines.
"""

from __future__ import annotations

import time
from abc import ABC, abstractmethod
from typing import Any

from ..context import RunContext
from ..manifest import StepRecord


class Step(ABC):
    name: str = ""
    inputs: tuple[str, ...] = ()
    outputs: tuple[str, ...] = ()

    @abstractmethod
    def run(self, ctx: RunContext, **inputs: Any) -> dict[str, Any]:
        """Do the work. Must return a dict whose keys match `outputs`."""

    def execute(self, ctx: RunContext, **inputs: Any) -> dict[str, Any]:
        """Wrapper: validate inputs, time the run, record to manifest."""
        missing = set(self.inputs) - set(inputs)
        if missing:
            raise ValueError(f"Step {self.name!r} missing inputs: {sorted(missing)}")

        started = time.time()
        result = self.run(ctx, **inputs)
        duration_ms = (time.time() - started) * 1000

        extra_keys = set(result) - set(self.outputs)
        missing_outs = set(self.outputs) - set(result)
        if extra_keys or missing_outs:
            raise ValueError(
                f"Step {self.name!r} bad outputs: "
                f"missing={sorted(missing_outs)} extra={sorted(extra_keys)}"
            )

        n_in = _row_count(inputs)
        n_out = _row_count(result)
        ctx.manifest.add(StepRecord(
            step=self.name,
            dataset=ctx.input.dataset,
            source=ctx.input.name,
            started_at=started,
            duration_ms=duration_ms,
            n_in=n_in,
            n_out=n_out,
        ))
        return result


def _row_count(d: dict) -> int | None:
    """Best-effort 'how many rows did this step see/produce'. Looks at
    the first list-like value; returns None if nothing list-like is present.
    Manifest readability nicety, not load-bearing."""
    for v in d.values():
        if isinstance(v, list):
            return len(v)
    return None

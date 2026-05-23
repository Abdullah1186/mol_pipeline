"""RunContext: per-run state passed to every step."""

from __future__ import annotations

from dataclasses import dataclass

from .config import Input, Params
from .manifest import Manifest


@dataclass
class RunContext:
    run_id: str
    params: Params
    manifest: Manifest
    input: Input   # the source DB currently being processed

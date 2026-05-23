"""Per-run manifest.

Every run writes one JSON file at artifacts/<run_id>/manifest.json
containing one record per (step, dataset, model). This is the lineage
record — what ran, with what params, on what, producing what.

Why a manifest matters: when a result looks off three weeks from now,
the manifest answers 'what was the input row count? was filter_odd_e on?
when did it run? how long did each step take?' without re-running anything.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path


@dataclass
class StepRecord:
    step: str
    dataset: str
    source: str                 # logical name from params.yaml inputs[].name
    started_at: float           # unix seconds
    duration_ms: float
    n_in: int | None            # input row count, if meaningful for this step
    n_out: int | None           # output row count, if meaningful
    extras: dict = field(default_factory=dict)  # step-specific small facts


@dataclass
class Manifest:
    run_id: str
    params_hash: str
    artifacts_dir: Path
    records: list[StepRecord] = field(default_factory=list)

    @property
    def path(self) -> Path:
        return self.artifacts_dir / self.run_id / "manifest.json"

    def add(self, record: StepRecord) -> None:
        self.records.append(record)
        self._flush()

    def _flush(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "run_id": self.run_id,
            "params_hash": self.params_hash,
            "records": [asdict(r) for r in self.records],
        }
        self.path.write_text(json.dumps(payload, indent=2))


def new_run_id() -> str:
    """Timestamped, human-sortable run id."""
    return time.strftime("%Y%m%d-%H%M%S")

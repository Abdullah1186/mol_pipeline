"""WriteMetrics: append one row to the metrics CSV.

Wraps filter.py:18-25 (header) and filter.py:135-140 (row). The runner
is responsible for creating the file + header before the first call;
this step only appends.

Columns must match the old filter.py exactly — the baseline parity gate
diffs against them.
"""

from __future__ import annotations

import csv
from typing import Any

from ..context import RunContext
from .base import Step


COLUMNS = [
    "Dataset", "Source", "Filters_applied",
    "Total", "Odd_e_before_filtering_VUN",
    "Valid", "Validity", "Unique", "Uniqueness", "Novel", "Novelty",
    "Odd_e_after_filtering_VUN", "Final_after_odd_e_filter",
]


class WriteMetrics(Step):
    name = "write_metrics"
    inputs = (
        "total", "odd_before",
        "valid_count", "validity", "unique_count", "uniqueness",
        "novel_count", "novelty",
        "odd_after", "final_kept", "filters_applied",
    )
    outputs = ("csv_path",)

    def run(self, ctx: RunContext, **inputs: Any) -> dict[str, Any]:
        path = ctx.params.metrics_csv
        with open(path, "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                ctx.input.dataset, ctx.input.name, inputs["filters_applied"],
                inputs["total"], inputs["odd_before"],
                inputs["valid_count"], inputs["validity"],
                inputs["unique_count"], inputs["uniqueness"],
                inputs["novel_count"], inputs["novelty"],
                inputs["odd_after"], inputs["final_kept"],
            ])
        return {"csv_path": path}


def write_header(path: str) -> None:
    """Create the CSV with header row. Called once by the runner before any step runs."""
    with open(path, "w", newline="") as f:
        csv.writer(f).writerow(COLUMNS)

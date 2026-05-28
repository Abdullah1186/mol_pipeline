"""Runner: build the step list and execute it per input.

What the runner does per input:
  1. load_db                              -> molecules, source_ids
  2. scan_odd_e_before(molecules)         -> odd_before
  3. evaluate(molecules, source_ids)      -> V/U/N + is_valid + is_unique
  4. apply_filters(...)                   -> kept_molecules, kept_source_ids
  5. scan_odd_e_after(kept_molecules)     -> odd_after (in final kept set)
  6. write_db(kept_source_ids)            -> output_path
  7. write_metrics(...)                   -> csv_path

Artifacts flow by *name*: each step declares its `inputs` tuple, the
runner pulls those names out of the per-pair `bag` dict.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from .config import Params, load_params
from .context import RunContext
from .manifest import Manifest, new_run_id
from .steps.apply_filters import ApplyFilters
from .steps.evaluate import Evaluate
from .steps.load_db import LoadDB
from .steps.scan_odd_e import ScanOddE
from .steps.write_db import WriteDB
from .steps.write_metrics import WriteMetrics, write_header


def _run_pair(ctx: RunContext) -> None:
    print(f"[{ctx.input.name}] starting")
    bag: dict = {}

    # 1. load
    bag.update(LoadDB().execute(ctx))
    total = len(bag["molecules"])

    # 2. odd-e scan on raw
    odd_before = ScanOddE("before").execute(ctx, molecules=bag["molecules"])["odd_count"]

    # 3. V/U/N measurement + per-mol flags
    bag.update(Evaluate().execute(
        ctx, molecules=bag["molecules"], source_ids=bag["source_ids"]
    ))

    # 4. apply user-selected filter criteria
    bag.update(ApplyFilters().execute(
        ctx,
        molecules=bag["molecules"], source_ids=bag["source_ids"],
        is_valid=bag["is_valid"], is_unique=bag["is_unique"],
    ))

    # 5. odd-e scan on the FINAL kept set (post-filter)
    odd_after = ScanOddE("after").execute(
        ctx, molecules=bag["kept_molecules"]
    )["odd_count"]

    # 6. write filtered DB
    WriteDB().execute(ctx, kept_source_ids=bag["kept_source_ids"])

    # 7. append CSV row
    WriteMetrics().execute(
        ctx,
        total=total,
        odd_before=odd_before,
        valid_count=bag["valid_count"],
        validity=bag["validity"],
        unique_count=bag["unique_count"],
        uniqueness=bag["uniqueness"],
        novel_count=bag["novel_count"],
        novelty=bag["novelty"],
        odd_after=odd_after,
        final_kept=len(bag["kept_molecules"]),
        filters_applied=ctx.params.filters.applied_label(),
    )
    print(f"[{ctx.input.name}] done — kept {len(bag['kept_molecules'])} / {total}")


def run(params: Params) -> Manifest:
    """Execute the pipeline against an in-memory Params. Returns the
    Manifest so callers (CLI, Streamlit UI, tests) can inspect lineage."""
    run_id = new_run_id()
    manifest = Manifest(
        run_id=run_id,
        params_hash=params.hash(),
        artifacts_dir=Path(params.artifacts_dir),
    )
    print(f"run_id={run_id}  params_hash={params.hash()}  n_inputs={len(params.inputs)}")

    # CSV header written once per run (matches old filter.py:18-25 — 'w' mode wipes existing)
    write_header(params.metrics_csv)

    for inp in params.inputs:
        ctx = RunContext(run_id=run_id, params=params, manifest=manifest, input=inp)
        _run_pair(ctx)

    print(f"manifest: {manifest.path}")
    print(f"metrics:  {params.metrics_csv}")
    return manifest


def main(params_path: str = "pipelines/params.yaml") -> None:
    """CLI entry: load YAML, run."""
    run(load_params(params_path))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", default="pipelines/params.yaml")
    args = parser.parse_args()
    main(args.params)

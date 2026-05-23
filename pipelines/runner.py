"""Runner: build the step list and execute it per (dataset, model).

This is a hand-rolled mini-DAG: a hardcoded linear order, no
topological sort, no caching, no parallelism. The point is to make the
*shape* of a pipeline obvious before letting DVC (Stage 2) hide it
behind YAML.

What the runner does per (dataset, model):
  1. load_db                              -> molecules, source_ids
  2. scan_odd_e_before(molecules)         -> odd_before
  3. evaluate(molecules, source_ids)      -> V/U/N + kept_*
  4. scan_odd_e_after(kept_molecules)     -> odd_after
  5. filter_odd_e(kept_*)                 -> kept_* (possibly fewer)
  6. write_db(kept_source_ids)            -> output_path
  7. write_metrics(... all numbers ...)   -> csv_path

Artifacts flow by *name*: each step declares its `inputs` tuple, the
runner pulls those names out of the per-pair `bag` dict. A missing or
typo'd name fails loudly via Step.execute()'s validation.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from .config import Params, load_params
from .context import RunContext
from .manifest import Manifest, new_run_id
from .steps.evaluate import Evaluate
from .steps.filter_odd_e import FilterOddE
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

    # 3. V/U/N + keep set
    eval_out = Evaluate().execute(
        ctx, molecules=bag["molecules"], source_ids=bag["source_ids"]
    )
    bag.update(eval_out)

    # 4. odd-e scan on kept set
    odd_after = ScanOddE("after").execute(
        ctx, molecules=bag["kept_molecules"]
    )["odd_count"]

    # 5. optional filter
    bag.update(FilterOddE().execute(
        ctx,
        kept_molecules=bag["kept_molecules"],
        kept_source_ids=bag["kept_source_ids"],
    ))

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

# The `pipelines/` package — a hand-rolled data pipeline

This is a small, deliberately minimal data pipeline that replaces the old
single-file `filter.py`. It exists for two reasons: (1) so this project is
easier to extend, and (2) so you can learn the vocabulary every real
pipeline tool uses (Airflow, Prefect, Dagster, Kedro, DVC) by reading ~200
lines of plain Python instead of a framework's docs.

---

## 1. What is a pipeline, and why?

`filter.py` used to be one 150-line script. Loading databases, building
molecules, computing metrics, filtering, writing the DB, and writing the
CSV were all interleaved in a single nested loop. That worked, but:

- You couldn't re-run just one part — say, recompute metrics without
  reloading 20k molecules.
- Adding a new dataset meant editing the loop, not a config file.
- Swapping a sink (write to Parquet instead of CSV) meant surgery on the
  middle of the script.
- There was no record of what produced any given output.

A **pipeline** is the same work split into named **steps** with declared
inputs and outputs, executed by a **runner**. That split is the whole
trick. Once you have it, every other pipeline feature (caching, lineage,
parallelism, retries, DAG viz) is a small addition on top.

---

## 2. The vocabulary

Each term is defined first abstractly, then grounded in *this* repo.

- **Step** — one named unit of work with declared inputs and outputs.
  → [steps/base.py](steps/base.py) defines the `Step` ABC. Every file in
  [steps/](steps/) is one Step.
- **Artifact** — a Step's named output. Can be in-memory (a list of
  molecules) or on-disk (a `.db` file).
  → `Step.outputs` lists the artifact names this Step produces.
- **Catalog** — owns naming conventions for pipeline outputs. Input
  paths are declared by the user in `params.yaml` (one entry per source
  DB), so the catalog itself stays small.
  → [catalog.py](catalog.py). `filtered_db_name(name, filter_odd_e=...)`
  is the only thing here now.
- **Params** — configuration separated from code. Flipping a knob means
  editing YAML, not Python.
  → [params.yaml](params.yaml) + the `Params` dataclass in
  [config.py](config.py).
- **Runner / DAG** — the thing that decides *what runs in what order*.
  In big tools this is a real graph the framework topologically sorts; in
  ours it's a hardcoded linear list of `.execute()` calls.
  → [runner.py](runner.py).
- **Run / Run ID** — one invocation of the pipeline. Everything that run
  produced is identified by its id.
  → `Manifest.run_id`, a sortable timestamp like `20260523-221124`.
- **Manifest / Lineage** — the per-run record of what ran, with what
  params, on what inputs, producing what.
  → [manifest.py](manifest.py) writes
  `artifacts/<run_id>/manifest.json`.
- **Idempotency / Caching** — re-running with the same inputs and params
  should be a no-op (or a cache hit). **Not implemented here** —
  deliberately deferred to Stage 2 (DVC). Mentioned so you know the gap
  exists.

---

## 3. How this pipeline is wired

```
LoadDB ──► ScanOddE("before") ──► Evaluate ──► ScanOddE("after") ──► FilterOddE ──► WriteDB
                                                                                  └─► WriteMetrics
```

Step-by-step:

| # | Step | Inputs | Outputs |
|---|------|--------|---------|
| 1 | `load_db` | — | `molecules`, `source_ids` |
| 2 | `scan_odd_e_before` | `molecules` | `odd_count` |
| 3 | `evaluate` | `molecules`, `source_ids` | V/U/N floats, V/U/N counts, `kept_molecules`, `kept_source_ids` |
| 4 | `scan_odd_e_after` | `kept_molecules` | `odd_count` |
| 5 | `filter_odd_e` | `kept_molecules`, `kept_source_ids` | same names (smaller lists if flag is on) |
| 6 | `write_db` | `kept_source_ids` | `output_path` |
| 7 | `write_metrics` | total + every metric + count | `csv_path` |

Artifacts flow by **name**: the runner pulls each Step's declared
`inputs` out of a shared dict, the Step validates them, runs, returns a
dict matching its declared `outputs`. A typo in either list fails loudly
the moment that Step is called — see `Step.execute()` in
[steps/base.py](steps/base.py).

---

## 4. Running it

```bash
# from the repo root
conda run -n geoldm python -m pipelines.runner
# or, via the back-compat shim:
conda run -n geoldm python filter.py
```

Outputs after a run:

- `molecular_metrics.csv` — one row per input DB. Columns: `Dataset,
  Source, Total, Odd_e_before_filtering_VUN, Valid, Validity, Unique,
  Uniqueness, Novel, Novelty, Odd_e_after_filtering_VUN,
  Final_after_odd_e_filter`.
- `<name>_filtered[_odd_e].db` — one ASE DB per input, in `output_dir`
  from [params.yaml](params.yaml). `<name>` is the `name:` field of the
  input entry.
- `artifacts/<run_id>/manifest.json` — per-step lineage record.

Fast iteration: edit [params.yaml](params.yaml) and trim the `inputs:`
list to one entry. One DB runs in seconds; the full list takes longer
proportional to size.

---

## 5. Extending it

### Add a new source DB

Just add one entry to `inputs:` in [params.yaml](params.yaml):

```yaml
inputs:
  - { path: /abs/path/to/any.db, dataset: qm9, name: my_run_42 }
```

`name` controls both the `Source` column in the CSV and the output
filename (`my_run_42_filtered.db`). No Python edits.

### Add a new dataset family (not qm9 / drugs)

If your DB needs a different atom encoder than the two we know about,
add one line to `_ALIASES` in [datasets.py](datasets.py) mapping your
`dataset:` value to the name `get_dataset_info` expects. Then use that
value in `params.yaml`.

### Add a new metric

1. Create `pipelines/steps/my_metric.py` subclassing `Step`, declaring
   its `inputs` and `outputs`.
2. Add one `MyMetric().execute(ctx, ...)` call inside `_run_pair()` in
   [runner.py](runner.py), wired to the right upstream artifact name(s).
3. If you want it in the CSV, extend `COLUMNS` and `WriteMetrics.inputs`
   in [steps/write_metrics.py](steps/write_metrics.py).

### Add a new sink (e.g., Parquet)

1. Mirror [steps/write_db.py](steps/write_db.py) as
   `steps/write_parquet.py`. Inputs whatever you want to persist.
2. Append `WriteParquet().execute(...)` to `_run_pair()` in the runner.

---

## 6. Where this maps onto real tools

Once you've understood the pieces here, every popular pipeline tool is a
re-skin of the same vocabulary:

| Our code | Airflow | Prefect | Dagster | DVC |
|----------|---------|---------|---------|-----|
| `Step` subclass | `Operator` / `@task` | `@task` | `@op` | one `stages:` entry |
| `Step.inputs`/`outputs` | XCom / explicit deps | task return values | `Out`/`In` | `deps:` / `outs:` |
| `catalog.py` | (per-operator paths) | (caller passes paths) | `IOManager` | `outs:` paths in `dvc.yaml` |
| `params.yaml` | Variables / DAG params | Blocks / parameters | `Config` | `params.yaml` (native) |
| `runner.py` (linear list) | DAG (topo-sorted graph) | flow graph | job/asset graph | DAG from `deps`/`outs` |
| `manifest.json` | run history DB + UI | run history / API | run timeline UI | `dvc.lock` |
| (nothing — re-runs always) | (re-runs unless skipped) | caching layer | asset materialization | `dvc repro` — true caching |

When you open an Airflow project later and see `Operator`, `DAG`, `XCom`
— you already know what they mean.

---

## 7. What's deliberately missing

- **Caching / idempotency.** Every run re-executes every step. Fixing
  this requires writing intermediate artifacts to disk with stable
  hashes — that's Stage 2 with DVC.
- **DAG topological sort.** Our runner is a hardcoded linear list. Real
  tools build the order from `inputs`/`outputs` automatically. We
  could've, but the linear list keeps the topology visible.
- **Parallelism.** Inputs are processed sequentially. They're
  independent and could run in parallel; not worth the complexity at
  this scale.
- **Retries / failure handling.** A step failure aborts the run. Fine
  for local batch work.
- **A UI.** Read `manifest.json` with `cat`.

These are the features that justify the existence of Airflow / Prefect /
Dagster. Now you know what they're buying you.

# Molecular Analysis Pipeline

A small, self-contained toolkit for filtering generated-molecule databases
and comparing their composition across models. Built around a hand-rolled
Python data pipeline and a Streamlit dashboard.

## What it does

Given one or more ASE `.db` files of generated molecules, the pipeline:

- Scores each molecule for **Validity** (RDKit sanitisable), **Uniqueness**
  (first-SMILES-wins), and **Novelty** (deferred).
- Lets you keep any combination of Valid / Unique / Even-electron
  molecules.
- Writes a filtered `.db`, a per-DB metrics CSV, and two JSON summaries
  per DB — `*.smiles.json` and `*.ecomp.json` (per-molecule elemental
  composition + a running total).
- The Streamlit UI drops the whole thing in a browser: upload DBs, pick
  filters, plot composition, element distributions, atom counts and
  molecular weights across datasets.

The dashboard plots read from `*.ecomp.json`, so comparing raw DBs is
instant and doesn't need the full pipeline to have run.

## Repo layout

```
app.py                    Streamlit UI entry point
filter.py                 3-line back-compat shim; calls pipelines.runner.main()

pipelines/                The data pipeline (CLI-driven, library-usable)
  README.md               Vocabulary + extension guide — read this to learn the
                          Step / Catalog / Params / Manifest pattern
  params.yaml             Local CLI config (ignored by the deployed UI)
  config.py               Params + FilterFlags dataclasses
  catalog.py              Output filename rules
  context.py              RunContext
  manifest.py             Per-run lineage record (artifacts/<run_id>/manifest.json)
  runner.py               Builds + executes the step list per input
  quick.py                One-shot "load → write raw JSON" helper used by the UI
  datasets.py             qm9/drugs atom encoder + atomic-number helpers
  steps/
    base.py               Step ABC: declared inputs/outputs, timed execute()
    load_db.py            ASE .db → list[(positions, atom_types)]
    scan_odd_e.py         Count odd-electron molecules
    evaluate.py           SMILES + V/U measurement, per-mol flags, cached SMILES
    apply_filters.py      AND the user-selected criteria into one mask
    write_db.py           Write kept molecules to a new ASE .db
    write_jsons.py        Emit smiles.json + ecomp.json
    write_metrics.py      Append a CSV row

plots/                    Matplotlib renderers (library, no Streamlit coupling)
  ecomp_bar.py            Average atom proportion per molecule (grouped bar)
  element_dist.py         Per-element count histogram, optional KDE overlay
  atoms_weight.py         Atom-count and molecular-weight KDE curves

baselines/                Reference CSVs for the pipeline parity gate
artifacts/                Per-run manifests (and CLI outputs under artifacts/cli/)
.uploads/                 Streamlit session dirs (one per browser session)

rdkit_functions.py        Original molecule-building + metrics module (reused as-is)
bond_analyze.py           Bond-order rules used by build_molecule
configs/datasets_config.py  Dataset metadata (atom_decoder, atomic_nb, …)

requirements.txt          Python dependencies
Dockerfile                Container image for Railway / Render / local docker
plot_functions.py         Reference implementations the plots/ modules were ported from
```

## Running

### Streamlit dashboard (local)

```bash
conda run -n geoldm streamlit run app.py --server.maxUploadSize=10000
```

Opens at <http://localhost:8501>. Three tabs:

1. **Upload** — drop `.db` files, assign a dataset family (`qm9` or
   `drugs`) and a label per file. Raw `ecomp.json` is auto-generated on
   drop so the dashboard works immediately.
2. **Filter** — tick V/U/Even-electrons in the sidebar, run the pipeline,
   view the metrics CSV + per-step manifest.
3. **Dashboard** — pick plots in the sidebar, pick datasets to compare,
   figures render automatically.

The `--server.maxUploadSize=10000` flag lifts Streamlit's 200 MB default
upload cap to 10 GB.

### CLI pipeline

```bash
conda run -n geoldm python -m pipelines.runner
```

Reads [pipelines/params.yaml](pipelines/params.yaml) and writes the
filtered `.db`, CSV, JSONs, and a manifest under `artifacts/cli/`.
`params.yaml` is for your own batch runs — the deployed UI never reads
it.

### Docker (local)

```bash
docker build -t mol-pipeline .
docker run -p 8501:8501 mol-pipeline
```

Same URL as the native run; identical to what Railway/Render will serve.

### Deploying to Railway / Render

Both hosts detect [Dockerfile](Dockerfile) automatically. Create a new
service from this repo, generate a public domain, done. No env vars
needed — the Dockerfile uses `$PORT`.

**Caveats:** the filesystem is ephemeral (sessions reset on restart),
RAM is small on free tiers, and the YAML `inputs:` paths point at a
local machine — only the UI is useful on a deployed instance.

## Environment

The pipeline + UI need `torch`, `rdkit`, `ase`, `matplotlib`, `scipy`,
`pandas`, `streamlit`, `pyyaml` — see [requirements.txt](requirements.txt).

Local setup (conda):

```bash
conda create -n geoldm python=3.9 -c conda-forge rdkit
conda activate geoldm
pip install -r requirements.txt
```

## How the pipeline works

Pipeline details — the Step / Artifact / Catalog / Manifest vocabulary,
how to add a new DB, how to add a new metric or sink — live in
[pipelines/README.md](pipelines/README.md). Read that one if you want to
extend the pipeline rather than just use it.

The one-line summary: each step has declared `inputs`/`outputs`,
[pipelines/runner.py](pipelines/runner.py) wires them together, every
run leaves a `manifest.json` behind.

## Development notes

- **Parity gate.** [baselines/](baselines/) holds reference CSVs. After
  any pipeline change, run `python -m pipelines.runner` and
  `diff baselines/molecular_metrics.baseline_off.csv
  artifacts/cli/molecular_metrics.csv`. Expected: empty diff.
- **Cleaning CLI outputs.** `rm -rf artifacts/cli/`. Everything the CLI
  writes clusters there.
- **Cleaning UI uploads.** `rm -rf .uploads/`. Or click "Clear uploads"
  in the sidebar for just the current session.
- **Plot defaults.** Figure sizes live inside each `plots/*.py` module.
  Dashboard tiles stack one per row at full container width.
- **SMILES caching.** `evaluate.py` emits `db_smiles` so both JSON
  writers reuse the SMILES instead of recomputing. Big perf win —
  dropped a 4-DB filter from ~6 min to ~2 min.

## Credits

Pipeline refactor + Streamlit UI built with Claude Code. Chemistry
primitives (`rdkit_functions.py`, `bond_analyze.py`, dataset configs)
are reused unchanged from the original codebase. The plotting helpers
in `plots/` are ports of
[plot_functions.py](plot_functions.py).

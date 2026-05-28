"""Streamlit UI for the molecule pipeline + plots.

Run with: streamlit run app.py

Flow:
  1. Drop .db files.
  2. Set dataset + name per file.
  3. Pick plots to render; for each plot pick raw or filtered source.
  4. (Only if any plot needs filtered) tweak filter settings.
  5. Run. Filter runs lazily (only if anything asked for filtered) and is
     cached per (db_path, filter_odd_e).
  6. Each plot renders one overlaid figure across all selected DBs, with
     a PNG download button.
"""

from __future__ import annotations

import io
import json
import shutil
import tempfile
import uuid
from pathlib import Path
from typing import Iterable

UPLOADS_ROOT = Path(".uploads")  # repo-local, gitignored

import pandas as pd
import streamlit as st
from matplotlib.figure import Figure

from pipelines import catalog
from pipelines.config import FilterFlags, Input, Params
from pipelines.datasets import info_for
from pipelines.runner import run as run_pipeline
from pipelines.steps.load_db import LoadDB
from pipelines.context import RunContext
from pipelines.manifest import Manifest, new_run_id
from plots.registry import PLOTS


RAW = FilterFlags(valid=False, unique=False, even_electrons=False)  # sentinel for "no filter at all"


st.set_page_config(page_title="Molecule Pipeline", layout="wide")
st.title("Molecule Pipeline")


# ---------------------------------------------------------------------------
# Session state
# ---------------------------------------------------------------------------
if "rows" not in st.session_state:
    st.session_state.rows = {}                # filename -> {dataset, name, path}
if "session_dir" not in st.session_state:
    UPLOADS_ROOT.mkdir(exist_ok=True)
    st.session_state.session_id = uuid.uuid4().hex[:8]
    st.session_state.session_dir = UPLOADS_ROOT / st.session_state.session_id
    st.session_state.session_dir.mkdir(exist_ok=True)
if "mol_cache" not in st.session_state:
    # (db_path, filter_odd_e) -> list[(positions, atom_types)]
    # filter_odd_e=None means "raw, no filter".
    st.session_state.mol_cache = {}
if "last_figures" not in st.session_state:
    st.session_state.last_figures = []        # list[(plot_name, Figure)]
if "last_filter_run" not in st.session_state:
    st.session_state.last_filter_run = None   # dict: csv_path, manifest_path, output_dir


def _wipe_session() -> None:
    """Delete this session's upload dir and reset all caches."""
    shutil.rmtree(st.session_state.session_dir, ignore_errors=True)
    st.session_state.session_dir.mkdir(parents=True, exist_ok=True)
    st.session_state.rows = {}
    st.session_state.mol_cache = {}
    st.session_state.last_figures = []
    st.session_state.last_filter_run = None


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
with st.sidebar:
    st.markdown("### Session")
    st.caption(f"id: `{st.session_state.session_id}`")
    st.caption(f"dir: `{st.session_state.session_dir}`")
    if st.button("Clear uploads", use_container_width=True):
        _wipe_session()
        st.rerun()


# ---------------------------------------------------------------------------
# 1. Upload
# ---------------------------------------------------------------------------
st.subheader("1. Drop ASE .db files")
uploaded = st.file_uploader(
    "Drop one or more .db files",
    type=["db"],
    accept_multiple_files=True,
    label_visibility="collapsed",
)
if uploaded:
    for f in uploaded:
        if f.name in st.session_state.rows:
            continue
        dest = st.session_state.session_dir / f.name
        dest.write_bytes(f.getbuffer())
        st.session_state.rows[f.name] = {
            "dataset": "qm9",
            "name": Path(f.name).stem,
            "path": dest,
        }


# ---------------------------------------------------------------------------
# 2. Per-file config
# ---------------------------------------------------------------------------
st.subheader("2. Per-file configuration")
if not st.session_state.rows:
    st.info("Drop files above to configure them.")
else:
    h1, h2, h3, h4 = st.columns([3, 2, 3, 1])
    h1.markdown("**File**")
    h2.markdown("**Dataset**")
    h3.markdown("**Name**")
    h4.markdown("**Remove**")
    for fname in list(st.session_state.rows):
        row = st.session_state.rows[fname]
        c1, c2, c3, c4 = st.columns([3, 2, 3, 1])
        c1.text(fname)
        row["dataset"] = c2.selectbox(
            "dataset", ["qm9", "drugs"],
            index=["qm9", "drugs"].index(row["dataset"]),
            key=f"ds_{fname}", label_visibility="collapsed",
        )
        row["name"] = c3.text_input(
            "name", value=row["name"],
            key=f"nm_{fname}", label_visibility="collapsed",
        )
        if c4.button("✕", key=f"rm_{fname}"):
            row["path"].unlink(missing_ok=True)
            del st.session_state.rows[fname]
            if not st.session_state.rows:
                _wipe_session()
            st.rerun()


# ---------------------------------------------------------------------------
# Helpers: load raw mols / filtered mols, both cached in session_state.
# Defined before section 3 so the filter-section cache-seeding can call them.
# ---------------------------------------------------------------------------
def _load_raw(path: Path, dataset: str, name: str) -> list:
    """Load all molecules from a DB, cached by path. Cache key uses RAW sentinel."""
    key = (str(path), RAW)
    if key not in st.session_state.mol_cache:
        inp = Input(path=str(path), dataset=dataset, name=name)
        run_id = new_run_id()
        manifest = Manifest(
            run_id=run_id,
            params_hash="ui",
            artifacts_dir=Path(st.session_state.session_dir) / "artifacts",
        )
        params = Params(
            filters=RAW, inputs=(inp,),
            output_dir=str(st.session_state.session_dir),
            metrics_csv=str(st.session_state.session_dir / "metrics_unused.csv"),
            artifacts_dir=str(st.session_state.session_dir / "artifacts"),
        )
        ctx = RunContext(run_id=run_id, params=params, manifest=manifest, input=inp)
        out = LoadDB().execute(ctx)
        st.session_state.mol_cache[key] = out["molecules"]
    return st.session_state.mol_cache[key]


def _load_filtered(rows: Iterable[dict], filters: FilterFlags) -> dict[str, list]:
    """Run the pipeline once over all rows whose filtered set isn't cached.
    Returns {name: kept_molecules} for every row."""
    missing = [
        r for r in rows
        if (str(r["path"]), filters) not in st.session_state.mol_cache
    ]
    if missing:
        run_dir = Path(tempfile.mkdtemp(prefix="run_", dir=st.session_state.session_dir))
        params = Params(
            filters=filters,
            inputs=tuple(
                Input(path=str(r["path"]), dataset=r["dataset"], name=r["name"])
                for r in missing
            ),
            output_dir=str(run_dir),
            metrics_csv=str(run_dir / "molecular_metrics.csv"),
            artifacts_dir=str(run_dir / "artifacts"),
        )
        run_pipeline(params)
        for r in missing:
            out_path = run_dir / catalog.filtered_db_name(r["name"], filters=filters)
            mols = _load_raw(out_path, r["dataset"], r["name"])
            st.session_state.mol_cache[(str(r["path"]), filters)] = mols
    return {
        r["name"]: st.session_state.mol_cache[(str(r["path"]), filters)]
        for r in rows
    }


# ---------------------------------------------------------------------------
# 3. Filter criteria + run + view CSV
#
# The FilterFlags chosen here are shared with the plot section below — if a
# plot asks for 'filtered' data, it gets the same flags.
# ---------------------------------------------------------------------------
st.subheader("3. Filter criteria")
fcol1, fcol2, fcol3 = st.columns(3)
f_valid = fcol1.checkbox("Validity", value=True, help="Keep RDKit-sanitisable mols.")
f_unique = fcol2.checkbox(
    "Uniqueness",
    value=True,
    disabled=not f_valid,
    help="Keep first occurrence per SMILES. Requires Validity.",
)
if not f_valid:
    f_unique = False
f_even = fcol3.checkbox(
    "Even electrons",
    value=False,
    help="Drop molecules whose total atomic number is odd.",
)
filters_selected = FilterFlags(valid=f_valid, unique=f_unique, even_electrons=f_even)
st.caption(f"Filters_applied = `{filters_selected.applied_label()}`")

filter_run_disabled = not st.session_state.rows
if filter_run_disabled:
    st.caption("Drop at least one file above to enable.")
filter_run_clicked = st.button(
    "Run filter pipeline",
    type="primary",
    disabled=filter_run_disabled,
    key="filter_section_run",
)

if filter_run_clicked:
    run_dir = Path(tempfile.mkdtemp(prefix="filter_", dir=st.session_state.session_dir))
    params = Params(
        filters=filters_selected,
        inputs=tuple(
            Input(path=str(r["path"]), dataset=r["dataset"], name=r["name"])
            for r in st.session_state.rows.values()
        ),
        output_dir=str(run_dir),
        metrics_csv=str(run_dir / "molecular_metrics.csv"),
        artifacts_dir=str(run_dir / "artifacts"),
    )
    with st.spinner(f"Running on {len(params.inputs)} input(s)…"):
        manifest = run_pipeline(params)
    # Seed plot cache so the same flags don't trigger a re-run from the plot section.
    for r in st.session_state.rows.values():
        out_path = run_dir / catalog.filtered_db_name(r["name"], filters=filters_selected)
        if out_path.exists():
            st.session_state.mol_cache[(str(r["path"]), filters_selected)] = (
                _load_raw(out_path, r["dataset"], r["name"])
            )
    st.session_state.last_filter_run = {
        "csv_path": params.metrics_csv,
        "manifest_path": str(manifest.path),
        "output_dir": str(run_dir),
        "filters_label": filters_selected.applied_label(),
    }
    st.success("Filter done.")

if st.session_state.last_filter_run is not None:
    lr = st.session_state.last_filter_run
    df = pd.read_csv(lr["csv_path"])
    st.dataframe(df, use_container_width=True)
    with st.expander("Manifest (per-step lineage)"):
        m = json.loads(Path(lr["manifest_path"]).read_text())
        st.caption(
            f"run_id={m['run_id']}  params_hash={m['params_hash']}  "
            f"records={len(m['records'])}  filters={lr['filters_label']}"
        )
        st.dataframe(pd.DataFrame(m["records"]), use_container_width=True)


# ---------------------------------------------------------------------------
# 4. Pick plots + per-plot source
# ---------------------------------------------------------------------------
st.subheader("4. Plots to render")
selected_plots: dict[str, str] = {}   # plot_name -> "raw" | "filtered"
for plot_name in PLOTS:
    col_a, col_b = st.columns([3, 2])
    enabled = col_a.checkbox(plot_name, key=f"plot_{plot_name}")
    if enabled:
        source = col_b.radio(
            "source", ["raw", "filtered"],
            horizontal=True,
            key=f"src_{plot_name}",
            label_visibility="collapsed",
        )
        selected_plots[plot_name] = source


# ---------------------------------------------------------------------------
# 5. Run plots
#
# Plots that pick 'filtered' use the FilterFlags from section 3.
# ---------------------------------------------------------------------------
st.subheader("5. Render plots")
needs_filter = any(src == "filtered" for src in selected_plots.values())
if needs_filter:
    st.caption(f"Plots marked 'filtered' will use: `{filters_selected.applied_label()}`")
run_disabled = not st.session_state.rows or not selected_plots
if run_disabled:
    st.caption("Drop at least one file and pick at least one plot to enable.")
run_clicked = st.button("Render plots", type="primary", disabled=run_disabled)


# ---------------------------------------------------------------------------
# 6. Execute on click
# ---------------------------------------------------------------------------
if run_clicked:
    rows = list(st.session_state.rows.values())
    figures: list[tuple[str, Figure]] = []

    # Pre-compute the data each plot needs. Group by source so filtered runs once.
    sources_needed = set(selected_plots.values())
    raw_data: dict[str, list] = {}
    filtered_data: dict[str, list] = {}

    with st.spinner("Loading…"):
        if "raw" in sources_needed:
            raw_data = {r["name"]: _load_raw(r["path"], r["dataset"], r["name"])
                        for r in rows}
        if "filtered" in sources_needed:
            filtered_data = _load_filtered(rows, filters_selected)

    # All loaded DBs share a dataset family in practice (one selectbox per row,
    # but for the dataset_info argument the plots need *a* dataset_info; for
    # mixed-dataset selections, use the first row's. Plots that need exact
    # decoder mappings should handle this themselves later.).
    info = info_for(rows[0]["dataset"])

    with st.spinner("Rendering plots…"):
        for plot_name, source in selected_plots.items():
            data = raw_data if source == "raw" else filtered_data
            fig = PLOTS[plot_name](data, info)
            figures.append((f"{plot_name}  ({source})", fig))

    st.session_state.last_figures = figures


# ---------------------------------------------------------------------------
# 7. Results
# ---------------------------------------------------------------------------
if st.session_state.last_figures:
    st.subheader("Results")
    for label, fig in st.session_state.last_figures:
        st.markdown(f"**{label}**")
        st.pyplot(fig)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        st.download_button(
            "Download PNG",
            data=buf.getvalue(),
            file_name=f"{label.replace(' ', '_').replace('(', '').replace(')', '')}.png",
            mime="image/png",
            key=f"dl_{label}",
        )
        st.divider()

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

UPLOADS_ROOT = Path(".uploads")  # repo-local, gitignored
MAX_UPLOAD_BYTES = 100 * 1024 * 1024  # 100 MB per file (deployed instance is RAM-constrained)

import pandas as pd
import streamlit as st

from pipelines.config import FilterFlags, Input, Params
from pipelines.runner import run as run_pipeline
from plots.ecomp_bar import bar_average_proportion


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
if "last_filter_run" not in st.session_state:
    st.session_state.last_filter_run = None   # dict: csv_path, manifest_path, output_dir


def _wipe_session() -> None:
    """Delete this session's upload dir and reset all caches."""
    shutil.rmtree(st.session_state.session_dir, ignore_errors=True)
    st.session_state.session_dir.mkdir(parents=True, exist_ok=True)
    st.session_state.rows = {}
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
        if f.size > MAX_UPLOAD_BYTES:
            st.error(
                f"{f.name} is {f.size / 1024 / 1024:.0f} MB — over the "
                f"{MAX_UPLOAD_BYTES // 1024 // 1024} MB per-file limit."
            )
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
# 4. Bar plot of average atom proportion (driven by ecomp.json files
#    produced by the filter pipeline). Pick any subset of the JSONs the
#    last run produced — raw + filtered for each DB — to overlay.
# ---------------------------------------------------------------------------
st.subheader("4. Average atom proportion (bar)")
if st.session_state.last_filter_run is None:
    st.caption("Run the filter pipeline above first — the bar plot reads its ecomp.json outputs.")
else:
    run_dir = Path(st.session_state.last_filter_run["output_dir"])
    json_paths = sorted(run_dir.glob("*.ecomp.json"))
    if not json_paths:
        st.caption("No ecomp.json files found in the last run's output dir.")
    else:
        # Label = filename stem (strip .ecomp); easier to read than full path.
        options = {p.name.replace(".ecomp.json", ""): str(p) for p in json_paths}
        picked = st.multiselect(
            "Datasets to compare",
            options=list(options),
            default=list(options),
            key="ecomp_bar_picked",
        )
        if st.button("Render bar plot", key="render_bar"):
            series = [(label, options[label]) for label in picked]
            with st.spinner("Rendering…"):
                fig = bar_average_proportion(series)
            st.pyplot(fig)
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
            st.download_button(
                "Download PNG",
                data=buf.getvalue(),
                file_name="ecomp_bar.png",
                mime="image/png",
                key="dl_ecomp_bar",
            )



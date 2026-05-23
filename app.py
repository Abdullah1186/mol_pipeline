"""Streamlit UI for the molecule pipeline.

Run with: streamlit run app.py

The pipeline is imported as a library — this file only handles file
uploads, building a Params object, calling pipelines.runner.run(), and
rendering results.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

from pipelines.config import Input, Params
from pipelines.runner import run as run_pipeline


st.set_page_config(page_title="Molecule Pipeline", layout="wide")
st.title("Molecule Pipeline")

# ---------------------------------------------------------------------------
# Session state: persist the per-file config rows + last-run results across
# Streamlit reruns (every widget interaction triggers a script rerun).
# ---------------------------------------------------------------------------
if "rows" not in st.session_state:
    st.session_state.rows = {}  # filename -> {"dataset": str, "name": str, "path": Path}
if "last_run" not in st.session_state:
    st.session_state.last_run = None  # dict with keys: csv_path, manifest_path, output_dir

# Each Streamlit session gets one temp dir; uploaded files live here for the
# lifetime of the session. Reusing the dir avoids duplicating large .db files
# on every rerun.
if "session_dir" not in st.session_state:
    st.session_state.session_dir = Path(tempfile.mkdtemp(prefix="mol_pipeline_ui_"))


# ---------------------------------------------------------------------------
# Sidebar: pipeline knobs + Run button
# ---------------------------------------------------------------------------
with st.sidebar:
    st.header("Settings")
    filter_odd_e = st.checkbox(
        "filter_odd_e",
        value=False,
        help="Drop odd-electron molecules from the final kept set.",
    )
    st.divider()
    run_disabled = len(st.session_state.rows) == 0
    run_clicked = st.button(
        "Run pipeline",
        type="primary",
        disabled=run_disabled,
        use_container_width=True,
    )
    if run_disabled:
        st.caption("Drop at least one .db file to enable.")


# ---------------------------------------------------------------------------
# 1. File uploader
# ---------------------------------------------------------------------------
st.subheader("1. Drop ASE .db files")
uploaded = st.file_uploader(
    "Drag and drop one or more .db files here",
    type=["db"],
    accept_multiple_files=True,
    label_visibility="collapsed",
)

# Persist newly-uploaded files to session_dir; add a row for each.
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
    # Header row
    h1, h2, h3, h4 = st.columns([3, 2, 3, 1])
    h1.markdown("**File**")
    h2.markdown("**Dataset**")
    h3.markdown("**Name** *(used in CSV + output filename)*")
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
            st.rerun()


# ---------------------------------------------------------------------------
# 3. Run
# ---------------------------------------------------------------------------
if run_clicked:
    # Each run gets its own output dir under the session dir, so successive
    # runs don't clobber each other and CLI runs in the repo aren't touched.
    run_dir = Path(tempfile.mkdtemp(prefix="run_", dir=st.session_state.session_dir))
    params = Params(
        filter_odd_e=filter_odd_e,
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
    st.session_state.last_run = {
        "csv_path": params.metrics_csv,
        "manifest_path": str(manifest.path),
        "output_dir": str(run_dir),
    }
    st.success("Done.")


# ---------------------------------------------------------------------------
# 4. Results
# ---------------------------------------------------------------------------
st.subheader("3. Last run")
if st.session_state.last_run is None:
    st.info("No run yet.")
else:
    lr = st.session_state.last_run
    df = pd.read_csv(lr["csv_path"])
    st.dataframe(df, use_container_width=True)

    with st.expander("Manifest (per-step lineage)"):
        manifest = json.loads(Path(lr["manifest_path"]).read_text())
        st.caption(
            f"run_id={manifest['run_id']}  "
            f"params_hash={manifest['params_hash']}  "
            f"records={len(manifest['records'])}"
        )
        st.dataframe(pd.DataFrame(manifest["records"]), use_container_width=True)

    with st.expander("Output files on disk"):
        for p in sorted(Path(lr["output_dir"]).rglob("*")):
            if p.is_file():
                st.code(str(p))

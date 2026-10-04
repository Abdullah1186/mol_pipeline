"""Molecular Analysis Pipeline — Streamlit dashboard.

Run with: streamlit run app.py

Layout:
  Sidebar: session controls, filter criteria + run, plot toggles,
           shared dataset picker, per-plot settings.
  Main:    file uploader, per-file config, metrics CSV, 2-column plot grid.
"""

from __future__ import annotations

import io
import json
import shutil
import tempfile
import uuid
from pathlib import Path

UPLOADS_ROOT = Path(".uploads")  # repo-local, gitignored

import pandas as pd
import streamlit as st

from pipelines.config import FilterFlags, Input, Params
from pipelines.quick import generate_raw_jsons
from pipelines.runner import run as run_pipeline
from plots.atoms_weight import atom_count_distribution, weight_distribution
from plots.ecomp_bar import bar_average_proportion
from plots.element_dist import _available_elements, element_distribution


st.set_page_config(page_title="Molecular Analysis Pipeline", layout="wide")
st.title("Molecular Analysis Pipeline")


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
if "raw_jsons" not in st.session_state:
    # (uploaded_filename, dataset, name) -> {"smiles": Path, "ecomp": Path}
    # Rebuilt whenever the (dataset, name) combo for a row changes.
    st.session_state.raw_jsons = {}
if "last_filter_run" not in st.session_state:
    st.session_state.last_filter_run = None   # dict: csv_path, manifest_path, output_dir


def _wipe_session() -> None:
    """Delete this session's upload dir and reset all caches."""
    shutil.rmtree(st.session_state.session_dir, ignore_errors=True)
    st.session_state.session_dir.mkdir(parents=True, exist_ok=True)
    st.session_state.rows = {}
    st.session_state.raw_jsons = {}
    st.session_state.last_filter_run = None


# ---------------------------------------------------------------------------
# Auto-gen for any (fname, dataset, name) combo lacking raw JSONs.
# Runs before the sidebar so the dataset multiselect can see fresh labels.
# ---------------------------------------------------------------------------
if st.session_state.rows:
    _raw_dir = st.session_state.session_dir / "raw"
    _todo = [
        (fname, row) for fname, row in st.session_state.rows.items()
        if (fname, row["dataset"], row["name"]) not in st.session_state.raw_jsons
    ]
    if _todo:
        _ptext = st.empty()
        _pbar = st.progress(0)
        for _i, (_fname, _row) in enumerate(_todo):
            _combo = (_fname, _row["dataset"], _row["name"])
            for _stale in [k for k in st.session_state.raw_jsons if k[0] == _fname]:
                _old = st.session_state.raw_jsons.pop(_stale)
                if _old["smiles"] is not None:
                    _old["smiles"].unlink(missing_ok=True)
                _old["ecomp"].unlink(missing_ok=True)
            _ptext.caption(
                f"Preparing raw JSONs for **{_row['name']}** ({_i + 1}/{len(_todo)})…"
            )
            _smi, _eco = generate_raw_jsons(
                Input(path=str(_row["path"]), dataset=_row["dataset"], name=_row["name"]),
                _raw_dir,
            )
            st.session_state.raw_jsons[_combo] = {"smiles": _smi, "ecomp": _eco}
            _pbar.progress((_i + 1) / len(_todo))
        _ptext.empty()
        _pbar.empty()
        # Rerun so the dashboard + sidebar multiselect pick up the fresh
        # entries immediately — without this, the user has to toggle a
        # control to force a rerun.
        st.rerun()


# ---------------------------------------------------------------------------
# Shared across sidebar + main pane: {label: ecomp_json_path}.
# ---------------------------------------------------------------------------
ecomp_options: dict[str, str] = {}
for _combo, _paths in st.session_state.raw_jsons.items():
    ecomp_options[f"{_combo[2]}_raw"] = str(_paths["ecomp"])
if st.session_state.last_filter_run is not None:
    _rd = Path(st.session_state.last_filter_run["output_dir"])
    for _p in sorted(_rd.glob("*_filtered_*.ecomp.json")):
        ecomp_options[_p.name.replace(".ecomp.json", "")] = str(_p)

# True when at least one dropped file hasn't had its raw JSONs generated yet.
# Used to disable plot controls so the user can't tweak settings that would
# produce an empty / inconsistent dashboard mid-upload.
_uploads_pending = any(
    (fname, row["dataset"], row["name"]) not in st.session_state.raw_jsons
    for fname, row in st.session_state.rows.items()
)
# A plot control is useful once uploads are done AND we actually have options.
_plot_controls_disabled = _uploads_pending or not ecomp_options


# ---------------------------------------------------------------------------
# Sidebar: session, filter criteria, plot toggles, dataset picker, per-plot
# settings. All in-pane plots read these values.
# ---------------------------------------------------------------------------
with st.sidebar:
    st.markdown("### Session")
    st.caption(f"id: `{st.session_state.session_id}`")
    st.caption(f"dir: `{st.session_state.session_dir}`")
    if st.button("Clear uploads", use_container_width=True):
        _wipe_session()
        st.rerun()

    st.divider()
    st.markdown("### Filter criteria")
    f_valid = st.checkbox("Validity", value=True,
                          help="Keep RDKit-sanitisable mols.")
    f_unique = st.checkbox("Uniqueness", value=True, disabled=not f_valid,
                           help="Keep first SMILES. Requires Validity.")
    if not f_valid:
        f_unique = False
    f_even = st.checkbox("Even electrons", value=False,
                         help="Drop odd-electron molecules.")
    filters_selected = FilterFlags(valid=f_valid, unique=f_unique,
                                   even_electrons=f_even)
    st.caption(f"applied = `{filters_selected.applied_label()}`")

    st.divider()
    st.markdown("### Plots")
    if _uploads_pending:
        st.caption("Preparing uploads — plots unlock when ready.")
    elif not ecomp_options:
        st.caption("Drop a .db file to enable.")
    # All checkboxes start unticked and stay disabled until uploads are done.
    show_bar = st.checkbox(
        "Average atom proportion (bar)",
        value=False, disabled=_plot_controls_disabled, key="dash_show_bar",
    )
    show_elem = st.checkbox(
        "Element distribution",
        value=False, disabled=_plot_controls_disabled, key="dash_show_elem",
    )
    show_atoms_wt = st.checkbox(
        "Atom count & molecular weight",
        value=False, disabled=_plot_controls_disabled, key="dash_show_atoms_wt",
    )

    st.markdown("### Datasets to compare")
    if _plot_controls_disabled:
        picked_labels: list[str] = []
        # Show a disabled, empty-looking placeholder so the layout doesn't jump.
        st.multiselect(
            "Datasets", options=[], default=[],
            disabled=True, key="dash_datasets_placeholder",
            label_visibility="collapsed",
        )
    else:
        picked_labels = st.multiselect(
            "Datasets",
            options=list(ecomp_options),
            default=list(ecomp_options),
            key="dash_datasets",
            label_visibility="collapsed",
        )

    # Per-plot settings (shown only when the relevant plot is enabled)
    picked_series = [(lab, ecomp_options[lab]) for lab in picked_labels]
    element_sym = None
    element_style = "histogram"
    if show_elem and not _plot_controls_disabled:
        st.markdown("### Element plot")
        _elements = _available_elements(picked_series) if picked_series else []
        if _elements:
            element_sym = st.selectbox("Element", _elements, key="dash_el_sym")
            element_style = st.radio(
                "Style", ["histogram", "histogram+kde"],
                horizontal=True, key="dash_el_style",
            )
        else:
            st.caption("Pick a dataset first.")

    atoms_wt_which = "Both"
    if show_atoms_wt and not _plot_controls_disabled:
        st.markdown("### Atom/weight plot")
        atoms_wt_which = st.radio(
            "Quantity", ["Atom count", "Molecular weight", "Both"],
            key="dash_aw_which",
        )



# ---------------------------------------------------------------------------
# Tabs: Filter | Dashboard. Both share st.session_state set up above.
# ---------------------------------------------------------------------------
tab_upload, tab_filter, tab_dashboard = st.tabs(["Upload", "Filter", "Dashboard"])

with tab_upload:
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
                # Drop this file's raw JSON entries + files from disk.
                for key in [k for k in st.session_state.raw_jsons if k[0] == fname]:
                    paths = st.session_state.raw_jsons.pop(key)
                    if paths["smiles"] is not None:
                        paths["smiles"].unlink(missing_ok=True)
                    paths["ecomp"].unlink(missing_ok=True)
                del st.session_state.rows[fname]
                if not st.session_state.rows:
                    _wipe_session()
                st.rerun()



with tab_filter:
    # ---------------------------------------------------------------------------
    # 3. Run filter pipeline + view metrics CSV.
    # Filter flags come from the sidebar.
    # ---------------------------------------------------------------------------
    st.subheader("3. Run filter pipeline")
    filter_run_disabled = not st.session_state.rows
    if filter_run_disabled:
        st.caption("Drop at least one file in the **Upload** tab to enable.")
    else:
        st.caption(f"Will apply: `{filters_selected.applied_label()}` (change in sidebar).")
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
        _rtext = st.empty()
        _rbar = st.progress(0)

        def _on_progress(i, total, name):
            if name is not None:
                _rtext.caption(f"Filtering **{name}** ({i + 1}/{total})…")
            _rbar.progress(i / total if total else 0)

        manifest = run_pipeline(params, on_progress=_on_progress)
        _rtext.empty()
        _rbar.empty()
        st.session_state.last_filter_run = {
            "csv_path": params.metrics_csv,
            "manifest_path": str(manifest.path),
            "output_dir": str(run_dir),
            "filters_label": filters_selected.applied_label(),
        }
        st.success("Filter done.")
        st.rerun()  # pick up the new filtered JSONs into ecomp_options

    if st.session_state.last_filter_run is not None:
        lr = st.session_state.last_filter_run
        with st.expander("Metrics CSV", expanded=True):
            st.dataframe(pd.read_csv(lr["csv_path"]), use_container_width=True)
        with st.expander("Manifest (per-step lineage)"):
            m = json.loads(Path(lr["manifest_path"]).read_text())
            st.caption(
                f"run_id={m['run_id']}  params_hash={m['params_hash']}  "
                f"records={len(m['records'])}  filters={lr['filters_label']}"
            )
            st.dataframe(pd.DataFrame(m["records"]), use_container_width=True)


    # ---------------------------------------------------------------------------
    # Dashboard: a 2-column grid of all plots the user enabled in the sidebar.
    # Plots auto-render whenever any sidebar control changes.
    # ---------------------------------------------------------------------------
    def _download_fig(fig, filename: str, key: str) -> None:
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        st.download_button(
            "Download PNG", data=buf.getvalue(),
            file_name=filename, mime="image/png", key=key,
        )


    def _render_bar() -> None:
        st.markdown("**Average atom proportion**")
        fig = bar_average_proportion(picked_series)
        st.pyplot(fig)
        _download_fig(fig, "ecomp_bar.png", "dl_bar")


    def _render_element() -> None:
        st.markdown(f"**Element distribution — {element_sym}**")
        fig = element_distribution(picked_series, element_sym, style=element_style)
        st.pyplot(fig)
        _download_fig(fig, f"element_dist_{element_sym}.png", "dl_el")


    def _render_atom_count() -> None:
        st.markdown("**Atom count per molecule**")
        fig = atom_count_distribution(picked_series)
        st.pyplot(fig)
        _download_fig(fig, "atom_count_kde.png", "dl_atoms")


    def _render_weight() -> None:
        st.markdown("**Molecular weight (g/mol)**")
        fig = weight_distribution(picked_series)
        st.pyplot(fig)
        _download_fig(fig, "molecular_weight_kde.png", "dl_weight")


    # Build the list of (title, render_fn) pairs the user has asked for.
    _tiles: list = []
    if show_bar and picked_series:
        _tiles.append(_render_bar)
    if show_elem and picked_series and element_sym is not None:
        _tiles.append(_render_element)
    if show_atoms_wt and picked_series:
        if atoms_wt_which in ("Atom count", "Both"):
            _tiles.append(_render_atom_count)
        if atoms_wt_which in ("Molecular weight", "Both"):
            _tiles.append(_render_weight)


with tab_dashboard:
    st.subheader("Dashboard")

    # Nudge the user over to the Filter tab if they haven't run it yet —
    # that's how `*_filtered_*` entries show up in the dataset multiselect.
    _has_filtered = any("_filtered_" in lab for lab in ecomp_options)
    if ecomp_options and not _has_filtered:
        st.info(
            "Run the filter pipeline in the **Filter** tab to see "
            "filtered datasets alongside the raw ones here."
        )

    if not ecomp_options:
        st.caption("Drop at least one .db file in the **Upload** tab to populate the dashboard.")
    elif not picked_series:
        st.caption("Pick at least one dataset in the sidebar.")
    elif not _tiles:
        st.caption("Enable at least one plot in the sidebar.")
    else:
        # 2-column grid, one tile per plot. Single tile goes full width.
        if len(_tiles) == 1:
            with st.container(border=True):
                _tiles[0]()
        else:
            for i in range(0, len(_tiles), 2):
                cols = st.columns(2)
                for col, render in zip(cols, _tiles[i:i + 2]):
                    with col, st.container(border=True):
                        render()



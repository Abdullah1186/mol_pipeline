"""One-shot helpers that skip the full pipeline.

Used by the UI to produce raw ecomp/smiles JSONs for a single dropped DB
without running Evaluate / ApplyFilters / WriteDB / WriteMetrics.

The full pipeline still exists — this is just a fast path for "give me
the raw JSONs so the bar plot has something to show before the user
touches Run".
"""

from __future__ import annotations

from pathlib import Path

from .config import FilterFlags, Input, Params
from .context import RunContext
from .manifest import Manifest, new_run_id
from .steps.load_db import LoadDB
from .steps.write_jsons import WriteJSONs


_NO_FILTERS = FilterFlags(valid=False, unique=False, even_electrons=False)


def generate_raw_jsons(inp: Input, out_dir: Path) -> tuple[Path, Path]:
    """Load `inp` and write <name>_raw.smiles.json + <name>_raw.ecomp.json
    into `out_dir`. Returns (smiles_path, ecomp_path)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    run_id = new_run_id()
    manifest = Manifest(
        run_id=run_id,
        params_hash="quick-raw",
        artifacts_dir=out_dir / "_quick_manifest",
    )
    params = Params(
        filters=_NO_FILTERS,
        inputs=(inp,),
        output_dir=str(out_dir),
        metrics_csv=str(out_dir / "_unused.csv"),
        artifacts_dir=str(out_dir / "_quick_manifest"),
    )
    ctx = RunContext(run_id=run_id, params=params, manifest=manifest, input=inp)

    load_out = LoadDB().execute(ctx)
    json_out = WriteJSONs("raw", basename=f"{inp.name}_raw").execute(
        ctx, molecules=load_out["molecules"]
    )
    return Path(json_out["smiles_path"]), Path(json_out["ecomp_path"])

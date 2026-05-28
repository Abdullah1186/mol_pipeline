"""Typed view of params.yaml.

Loading config is split from using it: the rest of the pipeline takes a
Params object, never reads the YAML file directly.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, asdict
from pathlib import Path

import yaml


@dataclass(frozen=True)
class FilterFlags:
    """Which criteria to AND together when deciding which molecules survive.

    - valid:          keep only molecules whose RDKit-sanitised SMILES is non-None
    - unique:         among valid mols, keep the first occurrence of each SMILES
                      (requires valid=True; meaningless otherwise)
    - even_electrons: drop molecules whose sum of atomic numbers is odd

    Novelty is intentionally absent — the pipeline doesn't load a
    reference SMILES set yet, so novelty is measured-as-0 and cannot be
    filtered on.
    """
    valid: bool = True
    unique: bool = True
    even_electrons: bool = False

    def applied_label(self) -> str:
        """Comma-joined names of the active filters, for the CSV row."""
        parts = []
        if self.valid: parts.append("valid")
        if self.unique: parts.append("unique")
        if self.even_electrons: parts.append("even_e")
        return ",".join(parts) if parts else "none"


@dataclass(frozen=True)
class Input:
    """One source DB the pipeline should process.

    - path:    absolute path to the source ASE .db
    - dataset: which dataset_info to use ('qm9' or 'drugs') — picks the
               atom encoder and odd-e atomic-number table
    - name:    label used in the CSV row and in output filenames
    """
    path: str
    dataset: str
    name: str


@dataclass(frozen=True)
class Params:
    filters: FilterFlags
    inputs: tuple[Input, ...]
    output_dir: str
    metrics_csv: str
    artifacts_dir: str

    def hash(self) -> str:
        blob = json.dumps(asdict(self), sort_keys=True).encode()
        return hashlib.sha256(blob).hexdigest()[:12]


def load_params(path: str | Path) -> Params:
    with open(path) as f:
        raw = yaml.safe_load(f)
    inputs = tuple(
        Input(path=str(i["path"]), dataset=str(i["dataset"]), name=str(i["name"]))
        for i in raw["inputs"]
    )
    f = raw.get("filters", {})
    return Params(
        filters=FilterFlags(
            valid=bool(f.get("valid", True)),
            unique=bool(f.get("unique", True)),
            even_electrons=bool(f.get("even_electrons", False)),
        ),
        inputs=inputs,
        output_dir=str(raw["output_dir"]),
        metrics_csv=str(raw["metrics_csv"]),
        artifacts_dir=str(raw["artifacts_dir"]),
    )

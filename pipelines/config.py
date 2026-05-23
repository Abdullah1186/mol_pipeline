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
    filter_odd_e: bool
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
    return Params(
        filter_odd_e=bool(raw["filter_odd_e"]),
        inputs=inputs,
        output_dir=str(raw["output_dir"]),
        metrics_csv=str(raw["metrics_csv"]),
        artifacts_dir=str(raw["artifacts_dir"]),
    )

"""Dataset metadata lookup.

filter.py mapped 'qm9' -> get_dataset_info('qm9'), everything else ->
get_dataset_info('geom'). Centralised here so steps don't repeat that
branch. Cached because get_dataset_info builds a non-trivial dict.
"""

from __future__ import annotations

from functools import lru_cache

from rdkit import Chem

from configs.datasets_config import get_dataset_info


_ALIASES = {"qm9": "qm9", "drugs": "geom"}


@lru_cache(maxsize=None)
def info_for(dataset: str) -> dict:
    try:
        canonical = _ALIASES[dataset]
    except KeyError:
        raise KeyError(f"Unknown dataset {dataset!r}. Known: {sorted(_ALIASES)}")
    return get_dataset_info(canonical, remove_h=False)


def atom_encoder_for(dataset: str) -> dict[str, int]:
    return info_for(dataset)["atom_encoder"]


def atomic_num_of(atom_idx: int, dataset_info: dict) -> int:
    """Atomic number for a dataset-local atom index. Lifted from filter.py:12."""
    if "atomic_nb" in dataset_info:
        return dataset_info["atomic_nb"][atom_idx]
    return Chem.GetPeriodicTable().GetAtomicNumber(dataset_info["atom_decoder"][atom_idx])


def has_odd_electrons(atom_types, dataset_info: dict) -> bool:
    """True iff the sum of atomic numbers is odd. Wraps the parity check
    repeated 3x in filter.py (lines 82-85, 119-122, 125-131)."""
    return sum(atomic_num_of(int(a), dataset_info) for a in atom_types) % 2 != 0

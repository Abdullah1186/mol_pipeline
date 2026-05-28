"""Output filename rules."""

from __future__ import annotations

from .config import FilterFlags


def filtered_db_name(name: str, *, filters: FilterFlags) -> str:
    """Filename for the filtered output DB.

    Suffix encodes which filters were applied so different runs don't
    overwrite each other:
      valid+unique               -> _filtered_VU.db
      valid+unique+even_e        -> _filtered_VUE.db
      etc. ('none' -> _filtered_raw.db)
    """
    parts = []
    if filters.valid: parts.append("V")
    if filters.unique: parts.append("U")
    if filters.even_electrons: parts.append("E")
    suffix = "".join(parts) if parts else "raw"
    return f"{name}_filtered_{suffix}.db"

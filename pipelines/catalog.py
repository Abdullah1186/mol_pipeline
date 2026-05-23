"""Output filename rules.

There used to be a hardcoded input-path registry here. Inputs now come
from params.yaml (each entry specifies its own absolute path), so this
module only owns the output-naming convention.
"""

from __future__ import annotations


def filtered_db_name(name: str, *, filter_odd_e: bool) -> str:
    """Filename for the filtered output DB."""
    suffix = "_filtered_odd_e" if filter_odd_e else "_filtered"
    return f"{name}{suffix}.db"

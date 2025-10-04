"""Default CLI parameters for target acquisition pipelines."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping


@dataclass(frozen=True)
class ModeDefaults:
    """Default values shared across CLI modes."""

    column: str
    chunk_size: int
    timeout: float
    limit: int | None
    offset: int


TARGET_MODE_DEFAULTS: Mapping[str, ModeDefaults] = {
    "chembl": ModeDefaults(
        column="target_chembl_id",
        chunk_size=5,
        timeout=30.0,
        limit=None,
        offset=0,
    ),
    "uniprot": ModeDefaults(
        column="uniprot_id",
        chunk_size=100,
        timeout=30.0,
        limit=None,
        offset=0,
    ),
    "iuphar": ModeDefaults(
        column="uniprot_id",
        chunk_size=100,
        timeout=30.0,
        limit=None,
        offset=0,
    ),
    "all": ModeDefaults(
        column="target_chembl_id",
        chunk_size=5,
        timeout=30.0,
        limit=None,
        offset=0,
    ),
}
"""Default CLI parameters per mode."""

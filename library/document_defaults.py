"""Default values shared across document CLI and configuration."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class DocumentModeDefaults:
    """Container for default values associated with a document pipeline mode."""

    column: str
    batch_size: int | None = None
    chunk_size: int | None = None
    sleep: float | None = None
    workers: int | None = None
    timeout: float | None = None


PUBMED_DEFAULTS = DocumentModeDefaults(
    column="PMID",
    batch_size=100,
    sleep=5.0,
    workers=1,
    timeout=10.0,
)

CHEMBL_DEFAULTS = DocumentModeDefaults(
    column="document_chembl_id",
    chunk_size=20,
    timeout=90.0,
)

ALL_DEFAULTS = DocumentModeDefaults(
    column="document_chembl_id",
    batch_size=50,
    chunk_size=20,
    sleep=5.0,
    workers=1,
    timeout=90.0,
)


__all__ = [
    "DocumentModeDefaults",
    "PUBMED_DEFAULTS",
    "CHEMBL_DEFAULTS",
    "ALL_DEFAULTS",
]

"""Default values shared across document-related CLIs and configuration."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PubmedDefaults:
    """Defaults for the PubMed document pipeline."""

    column: str = "PMID"
    sleep: float = 5.0
    workers: int = 1
    batch_size: int = 100
    chunk_size: int = 100
    timeout: float = 10.0
    limit: int | None = None


@dataclass(frozen=True)
class ChemblDefaults:
    """Defaults for the ChEMBL document pipeline."""

    column: str = "document_chembl_id"
    chunk_size: int = 5
    timeout: float = 30.0
    batch_size: int = 5
    sleep: float = 0.0
    workers: int = 1
    limit: int | None = None


@dataclass(frozen=True)
class AllDefaults:
    """Defaults for the combined ChEMBL + PubMed pipeline."""

    column: str = "document_chembl_id"
    chunk_size: int = 5
    sleep: float = 5.0
    workers: int = 1
    batch_size: int = 50
    timeout: float = 30.0
    chembl_batch_size: int = 5
    chembl_sleep: float = 0.0
    chembl_workers: int = 1
    pubmed_chunk_size: int = 100
    pubmed_timeout: float = 10.0
    limit: int | None = None


PUBMED_DEFAULTS = PubmedDefaults()
CHEMBL_DEFAULTS = ChemblDefaults()
ALL_DEFAULTS = AllDefaults()


MODE_DEFAULTS = {
    "pubmed": PUBMED_DEFAULTS,
    "chembl": CHEMBL_DEFAULTS,
    "all": ALL_DEFAULTS,
}

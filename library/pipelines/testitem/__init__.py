"""Test item pipeline orchestration helpers.

Changelog
========
* Add compatibility re-exports for the modern test item pipeline API.
"""

from __future__ import annotations

# ===== Modules =====
from .enrichment import enrich

# ===== Compatibility Exports =====
from library.testitem_pipeline.pubchem import (
    PUBCHEM_CID_CACHE_ENCODING,
    PUBCHEM_COLUMNS,
)

__all__ = [
    "enrich",
    "PUBCHEM_CID_CACHE_ENCODING",
    "PUBCHEM_COLUMNS",
]

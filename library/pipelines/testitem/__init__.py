"""Test item pipeline orchestration helpers.

Changelog
========
* Add compatibility re-exports for the modern test item pipeline API.
* Provide graceful fallbacks when the optional modern pipeline is unavailable.
"""

from __future__ import annotations

# ===== Modules =====
from .enrichment import enrich

# ===== Compatibility Exports =====
try:  # pragma: no cover - exercised indirectly via import side effects
    from library.testitem_pipeline.pubchem import (
        PUBCHEM_CID_CACHE_ENCODING,
        PUBCHEM_COLUMNS,
    )
except ModuleNotFoundError:  # pragma: no cover - platform specific dependency gap
    from library.common.text_utils import UTF8_ENCODING
    from library.pipelines.assay.chembl_assay import TESTITEM_PUBCHEM_COLUMNS

    PUBCHEM_CID_CACHE_ENCODING = UTF8_ENCODING
    PUBCHEM_COLUMNS = TESTITEM_PUBCHEM_COLUMNS

__all__ = [
    "enrich",
    "PUBCHEM_CID_CACHE_ENCODING",
    "PUBCHEM_COLUMNS",
]

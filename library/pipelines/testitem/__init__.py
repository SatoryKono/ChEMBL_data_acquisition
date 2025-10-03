"""Test item pipeline orchestration helpers.

Changelog
========
* Add compatibility re-exports for the modern test item pipeline API.
"""

from __future__ import annotations

# ===== Modules =====
from importlib import import_module
from importlib.util import find_spec

from .enrichment import enrich
from library.pipelines.assay.chembl_assay import TESTITEM_PUBCHEM_COLUMNS


# ===== Compatibility Exports =====
_PUBCHEM_COMPAT_MODULE = "library.testitem_pipeline.pubchem"

if find_spec(_PUBCHEM_COMPAT_MODULE) is not None:
    pubchem_module = import_module(_PUBCHEM_COMPAT_MODULE)
    PUBCHEM_CID_CACHE_ENCODING = pubchem_module.PUBCHEM_CID_CACHE_ENCODING
    PUBCHEM_COLUMNS = list(pubchem_module.PUBCHEM_COLUMNS)
else:
    PUBCHEM_CID_CACHE_ENCODING = "utf-8"
    PUBCHEM_COLUMNS = list(TESTITEM_PUBCHEM_COLUMNS)

__all__ = [
    "enrich",
    "PUBCHEM_CID_CACHE_ENCODING",
    "PUBCHEM_COLUMNS",
]

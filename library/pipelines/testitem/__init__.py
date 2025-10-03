"""Test item pipeline orchestration helpers.

Changelog
========
* Add compatibility re-exports for the modern test item pipeline API.
* Re-export public helpers from :mod:`library.testitem_pipeline` for legacy
  imports.
"""

from __future__ import annotations

# ===== Modules =====
from importlib import import_module
from typing import Any

from .enrichment import enrich

_PIPELINE_EXPORTS = {
    "PARENT_LOOKUP_SOURCE_CACHE",
    "PARENT_LOOKUP_SOURCE_LOOKUP",
    "PARENT_LOOKUP_SOURCE_PARTIAL",
    "PARENT_LOOKUP_SOURCE_SKIPPED",
    "PARENT_LOOKUP_SOURCE_SYNC",
    "ReadInputIdsResult",
    "TestitemPipelineOptions",
    "_DEFAULT_CATALOG_CFG",
    "_FETCH_ERROR_SAMPLE_SIZE",
    "_MOLECULE_HIERARCHY_COLUMNS",
    "_PUBCHEM_CACHE_SCHEMA_VERSION",
    "_TYPO_PARENT_COLUMN",
    "_prepare_pubchem_api_cfg",
    "_write_pubchem_cid_cache",
    "analyze_table_quality",
    "ensure_no_parant_column",
    "file_sha256",
    "fetch_testitems",
    "integrate_missing_identifiers",
    "load_parent_catalog",
    "query_parent_catalog",
    "read_input_ids",
    "run_testitem_pipeline",
    "update_parent_catalog_cache",
    "write_meta_yaml",
    "write_parent_catalog_cache",
}

_PUBCHEM_EXPORTS = {
    "PUBCHEM_CID_CACHE_ENCODING": "PUBCHEM_CID_CACHE_ENCODING",
    "PUBCHEM_COLUMNS": "PUBCHEM_COLUMNS",
}

__all__ = [
    "enrich",
    "PUBCHEM_CID_CACHE_ENCODING",
    "PUBCHEM_COLUMNS",
    *sorted(_PIPELINE_EXPORTS),
]


def __getattr__(name: str) -> Any:
    """Proxy access to the modern test item pipeline helpers."""

    if name in _PIPELINE_EXPORTS:
        module = import_module("library.testitem_pipeline")
        value = getattr(module, name)
        globals()[name] = value
        return value
    if name in _PUBCHEM_EXPORTS:
        module = import_module("library.testitem_pipeline.pubchem")
        value = getattr(module, _PUBCHEM_EXPORTS[name])
        if name == "PUBCHEM_COLUMNS":
            value = list(value)
        globals()[name] = value
        return value
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


def __dir__() -> list[str]:
    """Expose lazily loaded pipeline exports for introspection."""

    return sorted({*globals(), *__all__, *_PIPELINE_EXPORTS, *_PUBCHEM_EXPORTS})

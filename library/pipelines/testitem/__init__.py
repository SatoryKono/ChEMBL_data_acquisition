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
from importlib.util import find_spec

from .enrichment import enrich
from library.pipelines.assay.chembl_assay import TESTITEM_PUBCHEM_COLUMNS

try:
    from library.testitem_pipeline import (
        PARENT_LOOKUP_SOURCE_CACHE,
        PARENT_LOOKUP_SOURCE_LOOKUP,
        PARENT_LOOKUP_SOURCE_PARTIAL,
        PARENT_LOOKUP_SOURCE_SKIPPED,
        PARENT_LOOKUP_SOURCE_SYNC,
        ReadInputIdsResult,
        TestitemPipelineOptions,
        _DEFAULT_CATALOG_CFG,
        _FETCH_ERROR_SAMPLE_SIZE,
        _MOLECULE_HIERARCHY_COLUMNS,
        _PUBCHEM_CACHE_SCHEMA_VERSION,
        _TYPO_PARENT_COLUMN,
        analyze_table_quality,
        ensure_no_parant_column,
        file_sha256,
        fetch_testitems,
        integrate_missing_identifiers,
        load_parent_catalog,
        query_parent_catalog,
        read_input_ids,
        run_testitem_pipeline,
        update_parent_catalog_cache,
        write_meta_yaml,
        write_parent_catalog_cache,
        _prepare_pubchem_api_cfg,
        _write_pubchem_cid_cache,
        PUBCHEM_CID_CACHE_ENCODING as _PIPELINE_PUBCHEM_CID_CACHE_ENCODING,
        PUBCHEM_COLUMNS as _PIPELINE_PUBCHEM_COLUMNS,
    )
except ModuleNotFoundError as exc:  # pragma: no cover - environment specific
    _PIPELINE_IMPORT_ERROR = exc
else:
    _PIPELINE_IMPORT_ERROR = None

if _PIPELINE_IMPORT_ERROR is not None:  # pragma: no cover - environment specific
    missing_detail = str(_PIPELINE_IMPORT_ERROR)
    msg = (
        "library.pipelines.testitem requires optional modules from "
        "library.testitem_pipeline. Ensure your installation includes "
        "'library.testitem_pipeline.cli' and related files. "
        f"Original error: {missing_detail}"
    )
    raise ModuleNotFoundError(msg) from _PIPELINE_IMPORT_ERROR


# ===== Compatibility Exports =====
_PUBCHEM_COMPAT_MODULE = "library.testitem_pipeline.pubchem"

try:
    if find_spec(_PUBCHEM_COMPAT_MODULE) is None:
        raise ModuleNotFoundError(_PUBCHEM_COMPAT_MODULE)
    pubchem_module = import_module(_PUBCHEM_COMPAT_MODULE)

    PUBCHEM_CID_CACHE_ENCODING = pubchem_module.PUBCHEM_CID_CACHE_ENCODING
    PUBCHEM_COLUMNS = list(pubchem_module.PUBCHEM_COLUMNS)

except ModuleNotFoundError:
    PUBCHEM_CID_CACHE_ENCODING = _PIPELINE_PUBCHEM_CID_CACHE_ENCODING
    PUBCHEM_COLUMNS = list(_PIPELINE_PUBCHEM_COLUMNS)
else:
    PUBCHEM_CID_CACHE_ENCODING = pubchem_module.PUBCHEM_CID_CACHE_ENCODING
    PUBCHEM_COLUMNS = list(pubchem_module.PUBCHEM_COLUMNS)

__all__ = [
    "enrich",
    "PUBCHEM_CID_CACHE_ENCODING",
    "PUBCHEM_COLUMNS",
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
]

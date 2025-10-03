"""Legacy compatibility facade for the test item pipeline.

Changelog
========
* Redirect all imports to :mod:`library.pipelines.testitem` and emit a
  deprecation warning to encourage the modern namespace.
"""

from __future__ import annotations

import json
from warnings import warn

import requests

from library.clients import pubchem as pc
from library.common.csv_utils import write_csv_chunks_deterministic as write_csv_deterministic
from library.common.log import logger
from library.common.metadata import file_sha256, write_meta_yaml
from library.integration import pubchem_library as pl
from library.integration.chembl_client import ChemblClient
from library.qa.table_quality import analyze_table_quality
from library.qa.validation import validate_testitems
from library.pipelines import testitem as _modern_pipeline

warn(
    "library.testitem_pipeline is deprecated; import from "
    "library.pipelines.testitem instead.",
    DeprecationWarning,
    stacklevel=2,
)

globals().update({name: getattr(_modern_pipeline, name) for name in _modern_pipeline.__all__})

testitem_enrichment = _modern_pipeline.testitem_enrichment

__all__ = list(_modern_pipeline.__all__)


def _load_local_module(module_name: str):  # pragma: no cover - compatibility shim
    """Provide legacy attribute expected by historical tests and callers."""

    qualified_name = f"{__name__}.{module_name}"
    msg = (
        f"{qualified_name} (expected {module_name}.py or {module_name}.pyc to be bundled)"
    )
    raise ModuleNotFoundError(msg)

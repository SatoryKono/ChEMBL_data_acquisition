"""Test item pipeline orchestration helpers.

Changelog
========
* Add compatibility re-exports for the modern test item pipeline API.
* Re-export public helpers from :mod:`library.testitem_pipeline` for legacy
  imports.
* Document :mod:`library.pipelines.testitem` as the canonical module while the
  legacy :mod:`library.testitem_pipeline` shim continues to proxy to this
  package.
"""

from __future__ import annotations

# ===== Modules =====
from pathlib import Path

from library.clients import pubchem as pc  # noqa: F401
from library.common.csv_utils import (
    write_csv_chunks_deterministic as write_csv_deterministic,  # noqa: F401
)
from library.common.log import logger  # noqa: F401
from library.config import Config
from library.integration import pubchem_library as pl  # noqa: F401
from library.integration.chembl_client import ChemblClient  # noqa: F401
from library.metadata import file_sha256, write_meta_yaml  # noqa: F401
from library.pipelines.assay.chembl_assay import TESTITEM_PUBCHEM_COLUMNS  # noqa: F401
from library.pipelines.common import PipelineRunResult
from library.table_quality import analyze_table_quality  # noqa: F401
from library.validation import validate_testitems  # noqa: F401

from . import enrichment as testitem_enrichment  # noqa: F401
from .cli import (
    RAW_INDEX_COLUMN,
    TestitemPipelineOptions,
    _normalise_output_labels,
    _extract_metadata_parameters,
    _write_primary_metadata,
    ensure_raw_index_column,
    stream_missing_placeholder_frames,
    run_testitem_pipeline,
)
from .core import *  # noqa: F401,F403
from .core import __all__ as _CORE_EXPORTS

__all__ = list(
    dict.fromkeys(
        [
            *_CORE_EXPORTS,
            "ChemblClient",
            "TESTITEM_PUBCHEM_COLUMNS",
            "TestitemPipelineOptions",
            "analyze_table_quality",
            "enrich",
            "ensure_raw_index_column",
            "file_sha256",
            "logger",
            "pc",
            "pl",
            "RAW_INDEX_COLUMN",
            "_normalise_output_labels",
            "_extract_metadata_parameters",
            "_write_primary_metadata",
            "run_pipeline",
            "run_testitem_pipeline",
            "stream_missing_placeholder_frames",
            "testitem_enrichment",
            "validate_testitems",
            "write_csv_deterministic",
            "write_meta_yaml",
        ]
    )
)


def run_pipeline(config: Config, options: TestitemPipelineOptions) -> PipelineRunResult:
    """Execute the test item pipeline and return a common result structure."""

    cfg = config.model_copy(deep=True)
    pipelines = cfg.sources.chembl.pipelines
    section = pipelines.testitem
    updates: dict[str, object] = {}
    if getattr(options, "limit", None) is not None:
        updates["limit"] = options.limit
    if getattr(options, "offset", None) is not None:
        updates["offset"] = options.offset
    if updates:
        pipelines.testitem = section.model_copy(update=updates)

    output_candidate = getattr(options, "output_csv", None)
    if output_candidate is not None:
        output_path = Path(output_candidate)
    else:
        output_path = Path(options.input_csv)

    exit_code, _ = run_testitem_pipeline(cfg, options)
    reason = None if exit_code == 0 else "pipeline_failed"
    written = None if exit_code != 0 else True
    return PipelineRunResult(
        exit_code=exit_code,
        output_path=output_path,
        executed=True,
        reason=reason,
        written=written,
    )

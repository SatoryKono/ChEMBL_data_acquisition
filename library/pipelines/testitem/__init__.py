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

from . import core as _core
from . import enrichment as testitem_enrichment  # noqa: F401
from .cli import TestitemPipelineOptions, run_testitem_pipeline

for _name in _core.__all__:
    globals()[_name] = getattr(_core, _name)
del _name

_shared_exports = list(dict.fromkeys(_core.__all__))

__all__ = list(
    dict.fromkeys(
        [
            *_shared_exports,
            "ChemblClient",
            "TESTITEM_PUBCHEM_COLUMNS",
            "analyze_table_quality",
            "enrich",
            "file_sha256",
            "logger",
            "molecule_catalog",
            "pc",
            "pl",
            "run_pipeline",
            "validate_testitems",
            "write_csv_deterministic",
            "write_meta_yaml",
        ]
    )
)

if "testitem_enrichment" not in __all__:
    __all__.append("testitem_enrichment")


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

    exit_code = run_testitem_pipeline(cfg, options)
    reason = None if exit_code == 0 else "pipeline_failed"
    written = None if exit_code != 0 else True
    return PipelineRunResult(
        exit_code=exit_code,
        output_path=output_path,
        executed=True,
        reason=reason,
        written=written,
    )

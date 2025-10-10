"""Utility libraries for data acquisition and processing.

This package exposes commonly used helper functions and submodules. The
original implementation used absolute imports which fail when the package is
executed as part of a larger project.  The imports are now explicitly relative
so that ``python`` can resolve them correctly regardless of the working
directory.
"""

from __future__ import annotations

from importlib import import_module
from typing import Any

try:
    from .common.csv_utils import (
        sha256_file,
        write_csv_chunks_deterministic,
        write_csv_deterministic,
    )
    from .common.logging_setup import Logger, LoggerConfig, configure_logger
    from .common.timing import log_duration
    from .config import Config, load_config
    from .parser_schema import CSVExportArgs
except ModuleNotFoundError as exc:  # pragma: no cover - import guard
    missing = exc.name or "dependency"
    raise ModuleNotFoundError(
        "library package requires additional dependencies. "
        f"Missing module: '{missing}'. Install project requirements, e.g.\n"
        "  pip install -r requirements.txt\n"
        "  # or\n"
        "  pip install -e .[dev]"
    ) from exc

_LAZY_SUBMODULES = {
    "io",
    "postprocess",
    "qa",
    "schemas",
    "validation",
}

_EXPORTS: dict[str, str] = {
    "parse_terms": "library.pipelines.document.type_terms:parse_terms",
    "REVIEW_TERMS": "library.pipelines.document.type_terms:REVIEW_TERMS",
    "EXPERIMENTAL_TERMS": "library.pipelines.document.type_terms:EXPERIMENTAL_TERMS",
    "UNKNOWN_TERMS": "library.pipelines.document.type_terms:UNKNOWN_TERMS",
    "compute_scores": "library.pipelines.document.type_classifier:compute_scores",
    "decide_label": "library.pipelines.document.type_classifier:decide_label",
    "organism_classification": "library.pipelines.target:organism_classification",
    "protein_classification": "library.pipelines.target:protein_classification",
    "testitem_enrichment": "library.pipelines.testitem.enrichment",
    "OrganismClassificationRules": "library.pipelines.target.organism_classification:OrganismClassificationRules",
    "add_cellularity": "library.pipelines.target.organism_classification:add_cellularity",
    "add_cellularity_smart": "library.pipelines.target.organism_classification:add_cellularity_smart",
    "classify_by_lineage": "library.pipelines.target.organism_classification:classify_by_lineage",
    "classify_record": "library.pipelines.target.organism_classification:classify_record",
    "normalize": "library.pipelines.target.organism_classification:normalize",
    "TYPE_MULTICELLULAR": "library.pipelines.target.organism_classification:TYPE_MULTICELLULAR",
    "TYPE_UNICELLULAR": "library.pipelines.target.organism_classification:TYPE_UNICELLULAR",
    "TYPE_VIRAL": "library.pipelines.target.organism_classification:TYPE_VIRAL",
    "SidecarErrors": "library.sidecar:SidecarErrors",
}


def __getattr__(name: str) -> Any:
    """Dynamically import exposed helpers and submodules on first access."""

    if name in _LAZY_SUBMODULES:
        module = import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    if name in _EXPORTS:
        target = _EXPORTS[name]
        module_name, _, attribute = target.partition(":")
        module = import_module(module_name)
        value = getattr(module, attribute) if attribute else module
        globals()[name] = value
        return value
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


__all__ = [
    "parse_terms",
    "REVIEW_TERMS",
    "EXPERIMENTAL_TERMS",
    "UNKNOWN_TERMS",
    "compute_scores",
    "decide_label",
    "postprocess",
    "io",
    "qa",
    "schemas",
    "validation",
    "Config",
    "load_config",
    "SidecarErrors",
    "write_csv_deterministic",
    "write_csv_chunks_deterministic",
    "sha256_file",
    "Logger",
    "LoggerConfig",
    "configure_logger",
    "CSVExportArgs",
    "log_duration",
    "organism_classification",
    "protein_classification",
    "testitem_enrichment",
    "OrganismClassificationRules",
    "add_cellularity",
    "add_cellularity_smart",
    "classify_by_lineage",
    "classify_record",
    "normalize",
    "TYPE_MULTICELLULAR",
    "TYPE_UNICELLULAR",
    "TYPE_VIRAL",
]


def __dir__() -> list[str]:
    """Return available attributes including lazily loaded submodules."""

    return sorted({*globals(), *__all__, *_LAZY_SUBMODULES})

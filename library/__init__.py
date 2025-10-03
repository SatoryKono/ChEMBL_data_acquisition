"""Utility libraries for data acquisition and processing.

The package re-exports frequently used helpers to provide a compact public
surface.  Importing some of the underlying modules is comparatively expensive
and can introduce circular dependencies when scripts import :mod:`library`
before bootstrapping :data:`sys.path`.  To keep the package importable in
minimal environments (for example when ``scripts/get_data.py`` adjusts
``sys.path`` at runtime) we lazily resolve all heavyweight modules the first
time their attributes are requested.
"""

from __future__ import annotations

from importlib import import_module
from typing import Any

from .common.csv_utils import (
    sha256_file,
    write_csv_chunks_deterministic,
    write_csv_deterministic,
)
from .common.logging_setup import Logger, LoggerConfig, configure_logger
from .common.timing import log_duration
from .config import Config, load_config
from .parser_schema import CSVExportArgs

_LAZY_SUBMODULES = {
    "io",
    "qa",
    "schemas",
    "validation",
}

_LAZY_EXPORTS: dict[str, tuple[str, str | None]] = {
    "compute_scores": ("library.pipelines.document.type_classifier", "compute_scores"),
    "decide_label": ("library.pipelines.document.type_classifier", "decide_label"),
    "parse_terms": ("library.pipelines.document.type_terms", "parse_terms"),
    "REVIEW_TERMS": ("library.pipelines.document.type_terms", "REVIEW_TERMS"),
    "EXPERIMENTAL_TERMS": ("library.pipelines.document.type_terms", "EXPERIMENTAL_TERMS"),
    "UNKNOWN_TERMS": ("library.pipelines.document.type_terms", "UNKNOWN_TERMS"),
    "organism_classification": ("library.pipelines.target", None),
    "protein_classification": ("library.pipelines.target", None),
    "testitem_enrichment": ("library.pipelines.testitem.enrichment", None),
    "OrganismClassificationRules": (
        "library.pipelines.target.organism_classification",
        "OrganismClassificationRules",
    ),
    "add_cellularity": (
        "library.pipelines.target.organism_classification",
        "add_cellularity",
    ),
    "add_cellularity_smart": (
        "library.pipelines.target.organism_classification",
        "add_cellularity_smart",
    ),
    "classify_by_lineage": (
        "library.pipelines.target.organism_classification",
        "classify_by_lineage",
    ),
    "classify_record": (
        "library.pipelines.target.organism_classification",
        "classify_record",
    ),
    "normalize": ("library.pipelines.target.organism_classification", "normalize"),
    "TYPE_MULTICELLULAR": (
        "library.pipelines.target.organism_classification",
        "TYPE_MULTICELLULAR",
    ),
    "TYPE_UNICELLULAR": (
        "library.pipelines.target.organism_classification",
        "TYPE_UNICELLULAR",
    ),
    "TYPE_VIRAL": ("library.pipelines.target.organism_classification", "TYPE_VIRAL"),
    "SidecarErrors": ("library.sidecar", "SidecarErrors"),
}


def _resolve_lazy_attribute(name: str) -> Any:
    module_name, attr_name = _LAZY_EXPORTS[name]
    module = import_module(module_name)
    value = module if attr_name is None else getattr(module, attr_name)
    globals()[name] = value
    return value


def __getattr__(name: str) -> Any:
    """Dynamically import heavy submodules or attributes on first access."""

    if name in _LAZY_SUBMODULES:
        module = import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    if name in _LAZY_EXPORTS:
        return _resolve_lazy_attribute(name)
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


__all__ = [
    "compute_scores",
    "decide_label",
    "parse_terms",
    "REVIEW_TERMS",
    "EXPERIMENTAL_TERMS",
    "UNKNOWN_TERMS",
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
    """Return available attributes including lazily loaded members."""

    return sorted({*globals(), *__all__, *_LAZY_SUBMODULES, *_LAZY_EXPORTS})

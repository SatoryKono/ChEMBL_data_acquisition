"""Target postprocessing pipeline."""

from __future__ import annotations

from pathlib import Path

from library.postprocessing.target import (
    set_default_search_dir as _set_target_search_dir,
)

from .schema import TARGET_SCHEMA, validate_targets
from .steps import (
    PIPELINE_STEPS,
    enrich_target_synonyms,
    finalize_target_records,
    normalize_target_fields,
    run_target_pipeline,
)

_DEFAULT_SEARCH_DIR: Path | None = None

__all__ = [
    "PIPELINE_STEPS",
    "TARGET_SCHEMA",
    "enrich_target_synonyms",
    "finalize_target_records",
    "normalize_target_fields",
    "run_target_pipeline",
    "validate_targets",
    "set_default_search_dir",
]


def set_default_search_dir(search_dir: Path | str | None) -> None:
    """Proxy helper aligning isoform defaults with pipeline configuration."""

    resolved = None if search_dir is None else Path(search_dir)
    globals()["_DEFAULT_SEARCH_DIR"] = resolved
    _set_target_search_dir(search_dir)

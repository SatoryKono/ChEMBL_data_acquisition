"""Target post-processing utilities."""

from __future__ import annotations

from pathlib import Path

from . import isoform as _isoform_module
from .export import (
    TARGETS_OBJECT_COLUMNS,
    TARGETS_OPTIONAL_COLUMNS,
    TARGETS_REQUIRED_COLUMNS,
    prepare_targets_for_schema,
)
from .isoform import *  # noqa: F401,F403 - re-export legacy isoform helpers
from .main import postprocess_target_table

_DEFAULT_SEARCH_DIR = _isoform_module._DEFAULT_SEARCH_DIR
_FALLBACK_SEARCH_DIR = Path(_DEFAULT_SEARCH_DIR)


def set_default_search_dir(search_dir: Path | str | None) -> None:
    """Override the default directory used to locate target exports."""

    resolved = _FALLBACK_SEARCH_DIR if search_dir is None else Path(search_dir)
    _isoform_module._DEFAULT_SEARCH_DIR = resolved
    globals()["_DEFAULT_SEARCH_DIR"] = resolved


__all__ = [
    *_isoform_module.__all__,
    "postprocess_target_table",
    "TARGETS_REQUIRED_COLUMNS",
    "TARGETS_OPTIONAL_COLUMNS",
    "TARGETS_OBJECT_COLUMNS",
    "prepare_targets_for_schema",
    "set_default_search_dir",
]

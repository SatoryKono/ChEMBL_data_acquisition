"""Target post-processing utilities."""

from __future__ import annotations

from . import isoform as _isoform_module
from .export import (
    TARGETS_OBJECT_COLUMNS,
    TARGETS_OPTIONAL_COLUMNS,
    TARGETS_REQUIRED_COLUMNS,
    prepare_targets_for_schema,
)
from .isoform import *  # noqa: F401,F403 - re-export legacy isoform helpers
from .main import postprocess_target_table

__all__ = [
    *_isoform_module.__all__,
    "postprocess_target_table",
    "TARGETS_REQUIRED_COLUMNS",
    "TARGETS_OPTIONAL_COLUMNS",
    "TARGETS_OBJECT_COLUMNS",
    "prepare_targets_for_schema",
]

_DEFAULT_SEARCH_DIR = _isoform_module._DEFAULT_SEARCH_DIR

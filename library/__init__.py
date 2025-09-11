"""Utility libraries for data acquisition and processing.

This package exposes commonly used helper functions and submodules. The
original implementation used absolute imports which fail when the package is
executed as part of a larger project.  The imports are now explicitly relative
so that ``python`` can resolve them correctly regardless of the working
directory.
"""

from . import io, validation
from .config import Config, load_config
from .csv_utils import sha256_file, write_csv_deterministic
from .document_type_classifier import compute_scores, decide_label
from .document_type_terms import (
    EXPERIMENTAL_TERMS,
    REVIEW_TERMS,
    UNKNOWN_TERMS,
    parse_terms,
)
from .logging_setup import Logger, LoggerConfig, configure_logger
from .sidecar import SidecarErrors

__all__ = [
    "parse_terms",
    "REVIEW_TERMS",
    "EXPERIMENTAL_TERMS",
    "UNKNOWN_TERMS",
    "compute_scores",
    "decide_label",
    "io",
    "validation",
    "Config",
    "load_config",
    "SidecarErrors",
    "write_csv_deterministic",
    "sha256_file",
    "Logger",
    "LoggerConfig",
    "configure_logger",
]

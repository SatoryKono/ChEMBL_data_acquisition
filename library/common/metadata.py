"""Compatibility re-export for metadata helpers."""

from __future__ import annotations

from .metadata_writer import (
    Stats,
    file_sha256,
    record_quality_failure,
    write_meta_yaml,
)

__all__ = [
    "Stats",
    "file_sha256",
    "record_quality_failure",
    "write_meta_yaml",
]

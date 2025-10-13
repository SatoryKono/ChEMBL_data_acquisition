"""Compatibility re-export for sidecar utilities."""

from __future__ import annotations

from .common.sidecar import SidecarErrors, resolve_failure_chunk_size

__all__ = ["SidecarErrors", "resolve_failure_chunk_size"]

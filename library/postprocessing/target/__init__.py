"""Organism-level target table post-processing translated from Power Query."""

from __future__ import annotations

from .main import postprocess_target_table

__all__ = ["postprocess_target_table"]

try:  # pragma: no cover - compatibility shim for legacy imports
    from ..isoform import process_targets as _process_targets
except Exception:  # pragma: no cover - defensive fallback
    pass
else:  # pragma: no cover - compatibility shim
    process_targets = _process_targets
    __all__.append("process_targets")

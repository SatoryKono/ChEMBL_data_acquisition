"""Utility helpers shared by post-processing step definitions."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from library.common.log import logger

__all__ = ["build_step_payload"]


def _normalize(value: Any) -> Any:
    """Recursively normalise ``value`` so it can be serialised safely."""

    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): _normalize(val) for key, val in value.items()}
    if isinstance(value, list):
        return [_normalize(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_normalize(item) for item in value)
    if isinstance(value, set):
        return sorted(_normalize(item) for item in value)
    return value


def build_step_payload(step: str, **params: Any) -> Mapping[str, Any]:
    """Return a lightweight description of a configured post-processing step."""

    normalised = {str(key): _normalize(val) for key, val in params.items()}
    logger.debug("postprocess_step_config", step=step, params=normalised)
    return {"step": step, "params": normalised}

"""Helpers for annotating exported tables with pipeline metadata."""

from __future__ import annotations

from datetime import UTC, datetime
from functools import lru_cache

import pandas as pd

from ...common.run_context import get_current
from ...project_version import get_pipeline_version as _get_pipeline_version

# Re-export ``get_pipeline_version`` for backwards compatibility with existing imports.
get_pipeline_version = _get_pipeline_version


@lru_cache(maxsize=1)
def get_timestamp_utc() -> str:
    """Return an ISO 8601 timestamp representing the pipeline execution time."""

    context = get_current()
    if context is not None and context.generated_at:
        return context.generated_at
    return datetime.now(UTC).isoformat()


@lru_cache(maxsize=1)
def pipeline_metadata() -> dict[str, str]:
    """Return a mapping with pipeline metadata columns."""

    return {
        "pipeline_version": get_pipeline_version(),
        "timestamp_utc": get_timestamp_utc(),
    }


def add_pipeline_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """Return ``df`` with pipeline metadata columns added."""

    if df.empty:
        # Preserve dtypes while ensuring the columns exist for empty frames.
        metadata_values = pipeline_metadata()
        result = df.copy()
        for column, _value in metadata_values.items():
            result[column] = pd.Series(dtype="string")
        return result

    result = df.copy()
    for column, value in pipeline_metadata().items():
        result[column] = value
    return result


__all__ = [
    "add_pipeline_metadata",
    "get_pipeline_version",
    "get_timestamp_utc",
    "pipeline_metadata",
]

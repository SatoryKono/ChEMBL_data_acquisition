"""Helpers for building correlation matrices from quality profiling."""

from __future__ import annotations

import warnings
from collections.abc import Sequence

import pandas as pd

from ..common.log import logger
from ..table_quality import TableQualityProfiler as _LegacyTableQualityProfiler
from .qc_report import (
    _TABLE_PROFILER_TYPES,
    TableQualityProfilerLike,
    _build_reports_from_profiler,
    _is_table_profiler_instance,
    _prepare_filtered_frame,
    _validate_table_name,
)

try:  # pragma: no cover - optional compatibility import
    from ..qa.table_quality import TableQualityProfiler as _QaTableQualityProfiler
except ImportError:  # pragma: no cover - qa profiler not available
    _QaTableQualityProfiler = None  # type: ignore[assignment]

if (
    _QaTableQualityProfiler is not None
    and _QaTableQualityProfiler is not _LegacyTableQualityProfiler
    and _QaTableQualityProfiler not in _TABLE_PROFILER_TYPES
):
    _TABLE_PROFILER_TYPES = (*_TABLE_PROFILER_TYPES, _QaTableQualityProfiler)

TableQualityProfiler = _LegacyTableQualityProfiler


def build_correlation_matrix(
    frame: pd.DataFrame | None,
    *,
    table_name: str = "table",
    include_columns: Sequence[str] | None = None,
    exclude_columns: Sequence[str] | None = None,
    sample_rows: int | None = None,
    method: str = "pearson",
    profiler: TableQualityProfilerLike | None = None,
) -> pd.DataFrame:
    """Return the numeric correlation matrix for ``frame`` without file output."""

    _validate_table_name(table_name)

    profiler_instance: TableQualityProfilerLike | None = None
    fallback_to_frame = False
    invalid_profiler_type: type | None = None

    if profiler is not None:
        if _is_table_profiler_instance(profiler):
            profiler_instance = profiler
        elif frame is None:
            raise TypeError("profiler must be a TableQualityProfiler instance")
        else:
            fallback_to_frame = True
            invalid_profiler_type = type(profiler)
            logger.warning(
                "table_quality_profiler_invalid",
                table_name=table_name,
                profiler_type=invalid_profiler_type.__name__,
            )

    if profiler_instance is None:
        if frame is None:
            raise ValueError("frame is required when profiler is not provided")

        filtered = _prepare_filtered_frame(
            frame,
            table_name=table_name,
            sample_rows=sample_rows,
            include_columns=include_columns,
            exclude_columns=exclude_columns,
        )

        profiler_instance = TableQualityProfiler()
        profiler_instance.consume(filtered)

        if fallback_to_frame:
            warnings.warn(
                "Ignoring incompatible profiler %r; falling back to frame data"
                % (
                    invalid_profiler_type.__name__
                    if invalid_profiler_type is not None
                    else "unknown"
                ),
                RuntimeWarning,
                stacklevel=2,
            )

    _, numeric_candidates = _build_reports_from_profiler(profiler_instance)
    if not numeric_candidates:
        return pd.DataFrame()

    return pd.DataFrame(numeric_candidates).corr(method=method)


def generate_correlation_report(
    frame: pd.DataFrame,
    *,
    table_name: str,
    include_columns: Sequence[str] | None = None,
    exclude_columns: Sequence[str] | None = None,
    sample_rows: int | None = None,
    method: str = "pearson",
    profiler: TableQualityProfilerLike | None = None,
) -> pd.DataFrame:
    """Return the numeric correlation matrix for ``frame`` using defaults."""

    return build_correlation_matrix(
        frame,
        table_name=table_name,
        include_columns=include_columns,
        exclude_columns=exclude_columns,
        sample_rows=sample_rows,
        method=method,
        profiler=profiler,
    )

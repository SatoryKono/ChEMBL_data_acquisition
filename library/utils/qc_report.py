"""Helpers for generating in-memory quality control summaries."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from ..common.log import logger
from ..table_quality import (
    TableQualityProfiler as _LegacyTableQualityProfiler,
)
from ..table_quality import (
    _apply_sampling_and_filters,
)

try:  # pragma: no cover - optional compatibility import
    from ..qa.table_quality import TableQualityProfiler as _QaTableQualityProfiler
except ImportError:  # pragma: no cover - qa table quality module absent
    _QaTableQualityProfiler = None  # type: ignore[assignment]

_TABLE_PROFILER_TYPES: tuple[type, ...]

if _QaTableQualityProfiler is not None and _QaTableQualityProfiler is not _LegacyTableQualityProfiler:
    _TABLE_PROFILER_TYPES = (_LegacyTableQualityProfiler, _QaTableQualityProfiler)
else:
    _TABLE_PROFILER_TYPES = (_LegacyTableQualityProfiler,)

if TYPE_CHECKING:  # pragma: no cover - type checker assistance only
    from ..qa.table_quality import TableQualityProfiler as _QaTableQualityProfilerType
else:  # pragma: no cover - runtime fallback for typing alias
    _QaTableQualityProfilerType = _LegacyTableQualityProfiler  # type: ignore[assignment]

TableQualityProfilerLike = _LegacyTableQualityProfiler | _QaTableQualityProfilerType

# Re-export the legacy profiler for typing/backwards compatibility.
TableQualityProfiler = _LegacyTableQualityProfiler


def _is_table_profiler_instance(candidate: object) -> bool:
    """Return ``True`` when ``candidate`` behaves like a table profiler."""

    if isinstance(candidate, _TABLE_PROFILER_TYPES):
        return True

    if candidate is None:
        return False

    required_attrs = ("_columns", "_accumulators")
    if not all(hasattr(candidate, attr) for attr in required_attrs):
        return False

    columns = candidate._columns
    accumulators = candidate._accumulators

    if not isinstance(columns, Sequence) or not isinstance(accumulators, Mapping):
        return False

    return True


def _validate_table_name(table_name: str) -> str:
    """Return the table label used for logging and validation."""
    table_label = Path(str(table_name)).name
    if not table_label:
        raise ValueError("table_name must contain at least one non-separator character")
    return table_label


def _prepare_filtered_frame(
    frame: pd.DataFrame,
    *,
    table_name: str,
    sample_rows: int | None,
    include_columns: Sequence[str] | None,
    exclude_columns: Sequence[str] | None,
) -> pd.DataFrame:
    """Apply sampling and column filters mirroring :mod:`table_quality`."""
    if not isinstance(frame, pd.DataFrame):
        raise TypeError("quality profiling expects a pandas DataFrame")

    include_tuple = tuple(include_columns) if include_columns is not None else None
    exclude_tuple = tuple(exclude_columns) if exclude_columns is not None else None

    filtered, _ = _apply_sampling_and_filters(
        frame,
        table_name=table_name,
        sample_rows=sample_rows,
        include_columns=include_tuple,
        exclude_columns=exclude_tuple,
        include_warning_logged=[False],
        exclude_warning_logged=[False],
        no_columns_logged=[False],
    )
    return filtered


def _build_reports_from_profiler(
    profiler: TableQualityProfilerLike,
) -> tuple[pd.DataFrame, dict[str, pd.Series]]:
    """Generate quality report rows and numeric candidates without file writes."""

    rows: list[dict[str, object]] = []
    numeric_candidates: dict[str, pd.Series] = {}
    for column in profiler._columns:  # pylint: disable=protected-access
        accumulator = profiler._accumulators[column]  # pylint: disable=protected-access
        row = accumulator.finalize()
        rows.append(row)
        if accumulator.numeric_cov >= 0.8:
            numeric_candidates[column] = pd.Series(accumulator.numeric_full, dtype=float)

    column_order = [
        "column",
        "non_null",
        "non_empty",
        "empty_pct",
        "unique_cnt",
        "unique_pct_of_non_empty",
        "pattern_cov_doi",
        "pattern_cov_issn",
        "pattern_cov_isbn",
        "pattern_cov_url",
        "pattern_cov_email",
        "bool_like_cov",
        "numeric_cov",
        "numeric_min",
        "numeric_p50",
        "numeric_p95",
        "numeric_max",
        "numeric_mean",
        "numeric_std",
        "date_cov",
        "date_min",
        "date_p50",
        "date_max",
        "text_len_min",
        "text_len_p50",
        "text_len_p95",
        "text_len_max",
        "guessed_roles",
        "top_values",
    ]
    quality_report = pd.DataFrame(rows, columns=column_order)
    return quality_report, numeric_candidates


def build_reports_from_profiler(
    profiler: TableQualityProfilerLike,
    *,
    correlation_method: str = "pearson",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return QC and correlation reports derived from ``profiler``.

    Parameters
    ----------
    profiler:
        Populated :class:`~library.table_quality.TableQualityProfiler` instance.
    correlation_method:
        Pandas correlation method forwarded to :meth:`pandas.DataFrame.corr`.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
        Pair of quality summary and correlation matrix DataFrames.
    """

    quality_report, numeric_candidates = _build_reports_from_profiler(profiler)
    if numeric_candidates:
        correlation_report = pd.DataFrame(numeric_candidates).corr(
            method=correlation_method
        )
    else:
        correlation_report = pd.DataFrame()
    return quality_report, correlation_report


def build_qc_summary(
    frame: pd.DataFrame | None,
    *,
    table_name: str = "table",
    include_columns: Sequence[str] | None = None,
    exclude_columns: Sequence[str] | None = None,
    sample_rows: int | None = None,
    profiler: TableQualityProfilerLike | None = None,
) -> pd.DataFrame:
    """Return the quality-control summary DataFrame for ``frame``."""

    _validate_table_name(table_name)

    profiler_instance: TableQualityProfilerLike | None = None
    if profiler is not None:
        if _is_table_profiler_instance(profiler):
            profiler_instance = profiler
        else:
            logger.warning(
                "table_quality_profiler_invalid",
                table_name=table_name,
                profiler_type=type(profiler).__name__,
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

    quality_report, _ = _build_reports_from_profiler(profiler_instance)
    return quality_report


def generate_qc_report(
    frame: pd.DataFrame,
    *,
    table_name: str,
    include_columns: Sequence[str] | None = None,
    exclude_columns: Sequence[str] | None = None,
    sample_rows: int | None = None,
    profiler: TableQualityProfilerLike | None = None,
) -> pd.DataFrame:
    """Return a deterministic quality-control report for ``frame``.

    The helper mirrors :func:`build_qc_summary` while providing a concise
    interface for callers that only require the default configuration.
    """

    return build_qc_summary(
        frame,
        table_name=table_name,
        include_columns=include_columns,
        exclude_columns=exclude_columns,
        sample_rows=sample_rows,
        profiler=profiler,
    )

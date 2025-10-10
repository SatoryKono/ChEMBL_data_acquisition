"""Helpers for generating in-memory quality control summaries."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Tuple

import pandas as pd

from ..table_quality import TableQualityProfiler, _apply_sampling_and_filters


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
    profiler: TableQualityProfiler,
) -> Tuple[pd.DataFrame, dict[str, pd.Series]]:
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


def build_qc_summary(
    frame: pd.DataFrame | None,
    *,
    table_name: str = "table",
    include_columns: Sequence[str] | None = None,
    exclude_columns: Sequence[str] | None = None,
    sample_rows: int | None = None,
    profiler: TableQualityProfiler | None = None,
) -> pd.DataFrame:
    """Return the quality-control summary DataFrame for ``frame``."""

    _validate_table_name(table_name)
    if profiler is not None:
        if not isinstance(profiler, TableQualityProfiler):
            raise TypeError("profiler must be a TableQualityProfiler instance")
        quality_report, _ = _build_reports_from_profiler(profiler)
        return quality_report

    if frame is None:
        raise ValueError("frame is required when profiler is not provided")

    filtered = _prepare_filtered_frame(
        frame,
        table_name=table_name,
        sample_rows=sample_rows,
        include_columns=include_columns,
        exclude_columns=exclude_columns,
    )

    profiler = TableQualityProfiler()
    profiler.consume(filtered)

    quality_report, _ = _build_reports_from_profiler(profiler)
    return quality_report

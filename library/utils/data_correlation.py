"""Helpers for building correlation matrices from quality profiling."""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from ..table_quality import TableQualityProfiler
from .qc_report import _build_reports_from_profiler, _prepare_filtered_frame, _validate_table_name


def build_correlation_matrix(
    frame: pd.DataFrame,
    *,
    table_name: str = "table",
    include_columns: Sequence[str] | None = None,
    exclude_columns: Sequence[str] | None = None,
    sample_rows: int | None = None,
    method: str = "pearson",
) -> pd.DataFrame:
    """Return the numeric correlation matrix for ``frame`` without file output."""

    _validate_table_name(table_name)
    filtered = _prepare_filtered_frame(
        frame,
        table_name=table_name,
        sample_rows=sample_rows,
        include_columns=include_columns,
        exclude_columns=exclude_columns,
    )

    profiler = TableQualityProfiler()
    profiler.consume(filtered)

    _, numeric_candidates = _build_reports_from_profiler(profiler)
    if not numeric_candidates:
        return pd.DataFrame()

    return pd.DataFrame(numeric_candidates).corr(method=method)

"""Helpers for the cell line data acquisition pipeline."""

from __future__ import annotations

from ...config import Config
from .chembl import (
    CELL_LINE_BASE_COLUMNS,
    CELL_LINE_COLUMN_ORDER,
    get_cell_lines,
)
from .pipeline import (
    CellLinePipelineOptions,
    CellLinePipelineResult,
    prepare_cellline_dataframe,
    read_cellline_identifiers,
    run_cellline_pipeline,
)

__all__ = [
    "CELL_LINE_BASE_COLUMNS",
    "CELL_LINE_COLUMN_ORDER",
    "CellLinePipelineOptions",
    "CellLinePipelineResult",
    "get_cell_lines",
    "prepare_cellline_dataframe",
    "read_cellline_identifiers",
    "run_cellline_pipeline",
    "run_pipeline",
]


def run_pipeline(
    config: Config, options: CellLinePipelineOptions
) -> CellLinePipelineResult:
    """Execute the cell line pipeline returning the detailed result."""

    return run_cellline_pipeline(config, options)

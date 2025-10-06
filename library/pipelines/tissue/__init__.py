"""Helpers for retrieving and processing ChEMBL tissue metadata."""

from ...config import Config
from .chembl import TISSUE_BASE_COLUMNS, TISSUE_COLUMN_ORDER, get_tissues
from .pipeline import (
    TissuePipelineOptions,
    TissuePipelineResult,
    prepare_tissue_dataframe,
    read_tissue_identifiers,
    run_tissue_pipeline,
)

__all__ = [
    "TISSUE_BASE_COLUMNS",
    "TISSUE_COLUMN_ORDER",
    "get_tissues",
    "prepare_tissue_dataframe",
    "read_tissue_identifiers",
    "run_tissue_pipeline",
    "run_pipeline",
    "TissuePipelineOptions",
    "TissuePipelineResult",
]


def run_pipeline(config: Config, options: TissuePipelineOptions) -> TissuePipelineResult:
    """Execute the tissue pipeline returning the detailed result."""

    return run_tissue_pipeline(config, options)

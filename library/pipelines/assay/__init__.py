"""ChEMBL assay pipeline helpers and orchestration utilities."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from ...config import Config
from .chembl_assay import (
    ACTIVITY_COLUMNS,
    ASSAY_COLUMNS,
    ASSAY_VARIANT_COLUMN_ALIASES,
    get_activities,
    get_assay,
    get_assays,
    get_testitem,
)
from .postprocessing import postprocess_assays, postprocess_file

__all__ = [
    "ACTIVITY_COLUMNS",
    "ASSAY_COLUMNS",
    "ASSAY_VARIANT_COLUMN_ALIASES",
    "AssayPipelineOptions",
    "get_activities",
    "get_assay",
    "get_assays",
    "get_testitem",
    "postprocess_assays",
    "postprocess_file",
    "run_pipeline",
]


@dataclass(slots=True)
class AssayPipelineOptions:
    """Simple configuration container for the assay pipeline."""

    input_csv: Path
    output_csv: Path
    limit: int | None = None
    offset: int = 0
    skip_existing: bool = False
    force: bool = False
    timeout: float | None = None
    batch_size: int | None = None


def run_pipeline(config: Config, options: AssayPipelineOptions) -> int:
    """Execute the assay pipeline using the CLI implementation."""

    from scripts import get_assay_data as _assay_cli

    args = argparse.Namespace(
        input_csv=Path(options.input_csv),
        output_csv=Path(options.output_csv),
        final_out=Path(options.output_csv),
        limit=options.limit,
        offset=options.offset,
        skip_existing=options.skip_existing,
        force=options.force,
        timeout=options.timeout,
        batch_size=options.batch_size,
    )
    return int(_assay_cli.run(config, args))

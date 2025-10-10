"""ChEMBL assay pipeline components and programmatic runner."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from library.config import Config
from library.pipelines.common import PipelineRunResult

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


@dataclass(slots=True)
class AssayPipelineOptions:
    """Parameters controlling programmatic assay pipeline execution.

    During orchestrated runs ``skip_existing`` mirrors the ``--skip-existing`` CLI
    flag so previously exported assay CSV files can be retained without rerunning
    the pipeline.
    """

    input_csv: Path
    output_csv: Path
    limit: int | None = None
    offset: int = 0
    timeout: float | None = None
    batch_size: int | None = None
    skip_existing: bool = False
    force: bool = False

def run_pipeline(config: Config, options: AssayPipelineOptions) -> PipelineRunResult:
    """Execute the assay pipeline with programmatic options."""

    from scripts import get_assay_data as assay_cli  # Lazy import to avoid cycles

    return assay_cli.run_assay_service(config, options)


__all__ = [
    "ACTIVITY_COLUMNS",
    "ASSAY_COLUMNS",
    "ASSAY_VARIANT_COLUMN_ALIASES",
    "AssayPipelineOptions",
    "PipelineRunResult",
    "get_activities",
    "get_assay",
    "get_assays",
    "get_testitem",
    "postprocess_assays",
    "postprocess_file",
    "run_pipeline",
]

"""ChEMBL assay pipeline components and programmatic runner."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

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
    """Parameters controlling programmatic assay pipeline execution."""

    input_csv: Path
    output_csv: Path
    limit: int | None = None
    offset: int = 0
    timeout: float | None = None
    batch_size: int | None = None
    skip_existing: bool = False
    force: bool = False


def _update_assay_config(
    cfg: Config,
    *,
    limit: int | None,
    offset: int,
    timeout: float | None,
    batch_size: int | None,
) -> None:
    pipelines = cfg.sources.chembl.pipelines
    section = pipelines.assay
    updates: dict[str, object] = {"offset": offset}
    if limit is not None:
        updates["limit"] = limit
    if timeout is not None:
        updates["timeout"] = timeout
    if batch_size is not None:
        updates["batch_size"] = batch_size
    pipelines.assay = section.model_copy(update=updates)


def run_pipeline(config: Config, options: AssayPipelineOptions) -> PipelineRunResult:
    """Execute the assay pipeline with programmatic options."""

    output_path = Path(options.output_csv)
    if options.skip_existing and output_path.exists() and not options.force:
        return PipelineRunResult(
            exit_code=0,
            output_path=output_path,
            executed=False,
            reason="skip_existing",
            written=False,
        )

    from scripts import get_assay_data as assay_cli  # Lazy import to avoid cycles

    cfg = config.model_copy(deep=True)
    _update_assay_config(
        cfg,
        limit=options.limit,
        offset=options.offset,
        timeout=options.timeout,
        batch_size=options.batch_size,
    )

    args = SimpleNamespace(
        input_csv=Path(options.input_csv),
        final_out=output_path,
        output_csv=output_path,
        limit=options.limit,
        offset=options.offset,
        timeout=options.timeout or cfg.assay.timeout,
        batch_size=options.batch_size or cfg.assay.batch_size,
        skip_existing=options.skip_existing,
        force=options.force,
    )

    exit_code = assay_cli.run(cfg, args)
    reason = None if exit_code == 0 else "pipeline_failed"
    written = None if exit_code != 0 else True
    return PipelineRunResult(
        exit_code=exit_code,
        output_path=output_path,
        executed=True,
        reason=reason,
        written=written,
    )


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

"""Helpers and orchestrators for the activity data pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from library.config import Config
from library.pipelines.activity.runner import (
    ActivityCommandOptions,
    MAX_ACTIVITY_CHUNK_SIZE,
    run_activity_pipeline,
)
from library.pipelines.common import PipelineRunResult

from .action_properties import (
    annotate_action_properties,
    build_activity_properties,
    infer_action_type,
)
from .activities import get_activities


@dataclass(slots=True)
class ActivityPipelineOptions:
    """Configuration forwarded to the activity CLI runner.

    The combined ``get_data`` command maps ``--skip-existing`` and ``--dry-run``
    to ``skip_existing`` and ``dry_run`` respectively so operators can either
    reuse existing exports or perform validation passes without writing output.
    """

    input_csv: Path
    output_csv: Path
    limit: int | None = None
    offset: int = 0
    timeout: float | None = None
    batch_size: int | None = None
    workers: int | None = None
    dry_run: bool = False
    skip_existing: bool = False
    force: bool = False


def _update_activity_config(
    cfg: Config,
    *,
    limit: int | None,
    offset: int,
    timeout: float | None,
    batch_size: int | None,
    workers: int | None,
) -> None:
    pipelines = cfg.sources.chembl.pipelines
    section = pipelines.activity
    updates: dict[str, object] = {"offset": offset}
    if limit is not None:
        updates["limit"] = limit
    if timeout is not None:
        updates["timeout"] = timeout
    if batch_size is not None:
        updates["batch_size"] = batch_size
    if workers is not None:
        updates["workers"] = workers
    pipelines.activity = section.model_copy(update=updates)


def run_pipeline(config: Config, options: ActivityPipelineOptions) -> PipelineRunResult:
    """Execute the activity pipeline using programmatic options."""

    output_path = Path(options.output_csv)
    if options.skip_existing and output_path.exists() and not options.force:
        return PipelineRunResult(
            exit_code=0,
            output_path=output_path,
            executed=False,
            reason="skip_existing",
            written=False,
        )

    cfg = config.model_copy(deep=True)
    _update_activity_config(
        cfg,
        limit=options.limit,
        offset=options.offset,
        timeout=options.timeout,
        batch_size=options.batch_size,
        workers=options.workers,
    )

    helper_options = ActivityCommandOptions(
        input_csv=Path(options.input_csv),
        output_csv=output_path,
        final_output=output_path,
        limit=options.limit,
        offset=options.offset,
        timeout=options.timeout,
        batch_size=options.batch_size,
        workers=options.workers,
        dry_run=options.dry_run,
        skip_existing=options.skip_existing,
        force=options.force,
    )

    exit_code = run_activity_pipeline(cfg, helper_options)
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
    "ActivityPipelineOptions",
    "PipelineRunResult",
    "annotate_action_properties",
    "build_activity_properties",
    "get_activities",
    "MAX_ACTIVITY_CHUNK_SIZE",
    "infer_action_type",
    "run_pipeline",
]

"""Helpers and orchestration shims for the activity data pipeline."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from ...config import Config
from .action_properties import (
    annotate_action_properties,
    build_activity_properties,
    infer_action_type,
)
from .activities import get_activities

__all__ = [
    "ActivityPipelineOptions",
    "annotate_action_properties",
    "build_activity_properties",
    "infer_action_type",
    "get_activities",
    "run_pipeline",
]


@dataclass(slots=True)
class ActivityPipelineOptions:
    """Configuration values required to execute the activity pipeline."""

    input_csv: Path
    output_csv: Path
    limit: int | None = None
    offset: int = 0
    skip_existing: bool = False
    force: bool = False
    dry_run: bool = False
    timeout: float | None = None
    batch_size: int | None = None
    workers: int | None = None


def run_pipeline(config: Config, options: ActivityPipelineOptions) -> int:
    """Execute the activity pipeline using the CLI implementation."""

    if options.dry_run:
        return 0

    from scripts import get_activity_data as _activity_cli

    args = argparse.Namespace(
        input_csv=Path(options.input_csv),
        output_csv=Path(options.output_csv),
        final_out=Path(options.output_csv),
        limit=options.limit,
        offset=options.offset,
        skip_existing=options.skip_existing,
        force=options.force,
        dry_run=options.dry_run,
        timeout=options.timeout,
        batch_size=options.batch_size,
        workers=options.workers,
    )
    return int(_activity_cli.run(config, args))

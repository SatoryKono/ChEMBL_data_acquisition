"""Helpers for executing the activity pipeline via reusable primitives."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd

from ...cli import Logger
from ...cli.pipeline_definition import PipelineDefinition
from ...cli_utils import run_pipeline
from ...common.fetch_retry import ChunkFailureTracker
from ...config import Config
from ..common import (
    ChunkedFetchConfig,
    CsvWriterConfig,
    PipelineRunResult,
    prepare_chunked_pipeline,
)

MetadataHook = Callable[[pd.DataFrame], pd.DataFrame]


def run_activity_pipeline(
    *,
    fetch_config: ChunkedFetchConfig,
    metadata_hooks: Sequence[MetadataHook],
    fetch_chunk: Callable[[Sequence[str]], pd.DataFrame],
    writer_config: CsvWriterConfig,
    definition_kwargs: Mapping[str, Any],
    cfg: Config,
    logger: Logger,
    output_path: Path,
    failure_path: Path,
    fetch_failure_path: Path,
    chunk_failures: ChunkFailureTracker,
) -> PipelineRunResult:
    """Execute the activity pipeline using shared chunked execution helpers."""

    fetcher, writer = prepare_chunked_pipeline(
        fetch_config=fetch_config,
        fetch_chunk=fetch_chunk,
        csv_writer=writer_config,
    )

    definition_params = dict(definition_kwargs)
    definition_params["metadata_hooks"] = tuple(metadata_hooks)
    definition_params["writer"] = writer
    definition = PipelineDefinition(**definition_params)

    try:
        exit_code = run_pipeline(
            definition=definition,
            fetcher=fetcher,
            output_path=output_path,
            failure_path=failure_path,
            cfg=cfg,
            logger=logger,
        )
    except Exception:
        logger.exception(
            "Activity pipeline execution failed during chunked processing."
        )
        chunk_failures.save(fetch_failure_path, cfg=cfg)
        raise
    else:
        chunk_failures.save(fetch_failure_path, cfg=cfg)

    reason = None if exit_code == 0 else "pipeline_failed"
    written = True if exit_code == 0 else None
    return PipelineRunResult(
        exit_code=exit_code,
        output_path=Path(output_path),
        executed=True,
        reason=reason,
        written=written,
    )


__all__ = ["run_activity_pipeline"]

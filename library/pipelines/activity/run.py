"""Helpers for executing the activity pipeline via reusable primitives."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd

from ...cli import Logger
from ...cli.pipeline_definition import PipelineDefinition
from ...cli_utils import PipelineExecutionResult, run_pipeline
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
    emit_legacy_artifacts: bool = True,
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

    execution: PipelineExecutionResult
    try:
        execution = run_pipeline(
            definition=definition,
            fetcher=fetcher,
            output_path=output_path,
            failure_path=failure_path,
            cfg=cfg,
            logger=logger,
            emit_legacy_artifacts=emit_legacy_artifacts,
        )
    except Exception:
        logger.exception(
            "Activity pipeline execution failed during chunked processing."
        )
        if emit_legacy_artifacts:
            chunk_failures.save(fetch_failure_path, cfg=cfg)
        else:
            fetch_failure_path.unlink(missing_ok=True)
            Path(f"{fetch_failure_path}.meta.yaml").unlink(missing_ok=True)
        raise
    else:
        if emit_legacy_artifacts:
            chunk_failures.save(fetch_failure_path, cfg=cfg)
        else:
            fetch_failure_path.unlink(missing_ok=True)
            Path(f"{fetch_failure_path}.meta.yaml").unlink(missing_ok=True)

    exit_code_attr = getattr(execution, "exit_code", None)
    exit_code = int(exit_code_attr if exit_code_attr is not None else execution)
    dataset_path = getattr(execution, "dataset_path", None) or output_path
    reason = None if exit_code == 0 else "pipeline_failed"
    written = True if exit_code == 0 else None
    return PipelineRunResult(
        exit_code=exit_code,
        output_path=Path(dataset_path),
        executed=True,
        reason=reason,
        written=written,
    )


__all__ = ["run_activity_pipeline"]

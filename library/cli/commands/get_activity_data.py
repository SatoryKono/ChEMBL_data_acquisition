from __future__ import annotations

import argparse
from collections.abc import Callable, Sequence
from dataclasses import dataclass
import sys
from pathlib import Path
from time import perf_counter

from library import io
from library.common.log import logger
from library.config import Config
from library.pipelines.assay.chembl_assay import MAX_ACTIVITY_CHUNK_SIZE

from . import _run

MIN_ACTIVITY_TIMEOUT = 60.0


@dataclass(slots=True)
class ActivityCommandOptions:
    """Options required to execute the activity pipeline programmatically."""

    input_csv: Path | str
    output_csv: Path | str | None = None
    final_output: Path | str | None = None
    limit: int | None = None
    offset: int = 0
    timeout: float | None = None
    batch_size: int | None = None
    workers: int | None = None
    dry_run: bool = False
    skip_existing: bool = False
    force: bool = False
    invocation: Sequence[str] | None = None


def _normalise_optional_path(value: Path | str | None) -> Path | None:
    """Return ``value`` as :class:`Path` when provided and meaningful."""

    if value in (None, argparse.SUPPRESS):
        return None
    return value if isinstance(value, Path) else Path(value)


def _normalise_path(value: Path | str) -> Path:
    """Convert a mandatory path-like value to :class:`Path`."""

    result = _normalise_optional_path(value)
    if result is None:  # pragma: no cover - defensive programming
        msg = "a concrete path value is required"
        raise ValueError(msg)
    return result


def run_activity_pipeline(
    cfg: Config,
    options: ActivityCommandOptions,
    *,
    runner: Callable[[Config, argparse.Namespace], int] | None = None,
    emit_completion_message: Callable[..., None] | None = None,
) -> int:
    """Execute the activity pipeline with orchestrator semantics.

    The helper mirrors :func:`scripts.get_activity_data.run` but accepts a
    structured ``options`` payload, enabling reuse from pipelines and tests.
    ``runner`` and ``emit_completion_message`` are primarily intended for
    dependency injection during testing; by default they fall back to the
    implementation from :mod:`scripts.get_activity_data`.
    """

    if runner is None or emit_completion_message is None:  # pragma: no cover - lazy import
        from scripts import get_activity_data as activity_cli

        if runner is None:
            runner = activity_cli.run_chembl
        if emit_completion_message is None:
            emit_completion_message = activity_cli._emit_completion_message

    if runner is None or emit_completion_message is None:  # pragma: no cover - defensive
        msg = "activity pipeline runner dependencies could not be resolved"
        raise RuntimeError(msg)

    start_time = perf_counter()

    active_logger = logger
    script_module = sys.modules.get("scripts.get_activity_data")
    if script_module is not None:
        patched_logger = getattr(script_module, "logger", None)
        if patched_logger is not None:
            active_logger = patched_logger

    args = argparse.Namespace(
        input_csv=_normalise_path(options.input_csv),
        skip_existing=bool(options.skip_existing),
        force=bool(options.force),
        dry_run=bool(options.dry_run),
    )
    args.limit = options.limit
    args.offset = options.offset
    args.timeout = (
        options.timeout
        if options.timeout is not None
        else getattr(getattr(cfg, "activity", object()), "timeout", None)
    )
    args.batch_size = (
        options.batch_size
        if options.batch_size is not None
        else getattr(getattr(cfg, "activity", object()), "batch_size", None)
    )
    args.workers = (
        options.workers
        if options.workers is not None
        else getattr(getattr(cfg, "activity", object()), "workers", None)
    )
    if options.invocation is not None:
        args.invocation = tuple(str(part) for part in options.invocation)

    final_output = _normalise_optional_path(options.final_output)
    legacy_output = _normalise_optional_path(options.output_csv)
    if final_output is None:
        if legacy_output is not None:
            output_path = legacy_output
        else:
            output_path = Path(io.default_output_path(args.input_csv, cfg.io))
        args.final_out = output_path
        args.output_csv = output_path
    else:
        output_path = final_output
        args.final_out = output_path
        args.output_csv = output_path

    batch_size = getattr(cfg.activity, "batch_size", None)
    if batch_size is not None and batch_size > MAX_ACTIVITY_CHUNK_SIZE:
        active_logger.warning(
            "activity_batch_size_clamped",
            configured=batch_size,
            limit=MAX_ACTIVITY_CHUNK_SIZE,
        )
        active_logger.warning(
            f"Configured batch size {batch_size} exceeds the hard cap of {MAX_ACTIVITY_CHUNK_SIZE}; "
            f"reducing to {MAX_ACTIVITY_CHUNK_SIZE}."
        )
        cfg.activity.batch_size = MAX_ACTIVITY_CHUNK_SIZE
        args.batch_size = MAX_ACTIVITY_CHUNK_SIZE

    timeout = getattr(cfg.activity, "timeout", None)
    if timeout is not None and timeout < MIN_ACTIVITY_TIMEOUT:
        active_logger.warning(
            "activity_timeout_clamped",
            configured=timeout,
            minimum=MIN_ACTIVITY_TIMEOUT,
        )
        active_logger.warning(
            f"Configured timeout {timeout} is below the minimum of {MIN_ACTIVITY_TIMEOUT}; "
            f"increasing to {MIN_ACTIVITY_TIMEOUT}."
        )
        minimum_timeout = float(MIN_ACTIVITY_TIMEOUT)
        cfg.activity.timeout = minimum_timeout
        try:
            args.timeout = float(minimum_timeout)
        except (TypeError, ValueError):  # pragma: no cover - defensive
            args.timeout = minimum_timeout

    retry_attempts = getattr(cfg.retry, "max_attempts", None)
    if retry_attempts is not None and retry_attempts <= 1:
        active_logger.warning(
            "activity_retry_disabled",
            configured=retry_attempts,
        )
        active_logger.warning(
            "Configured system.retry.max_attempts=%s disables urllib3 retry handling; "
            "increase to at least 2 to tolerate transient network issues.",
            retry_attempts,
        )

    api_retries = getattr(cfg.api, "retries", None)
    if api_retries is not None and api_retries <= 0:
        active_logger.warning(
            "activity_api_retry_disabled",
            configured=api_retries,
        )
        active_logger.warning(
            "Configured chembl.api.retries=%s disables client-level request retries; "
            "increase to 1 or more for resilience.",
            api_retries,
        )

    if args.skip_existing and output_path.exists() and not args.force:
        active_logger.info("pipeline_skip_existing", output=str(output_path))
        active_logger.info(
            f"Skipping execution because '{output_path}' already exists and --force was not provided."
        )
        emit_completion_message(
            output_path=output_path,
            processed_rows=None,
            duration_s=perf_counter() - start_time,
            mode="skip_existing",
        )
        return 0

    return runner(cfg, args)


def main(argv: Sequence[str] | None = None) -> int:
    """Run scripts.get_activity_data.

    Parameters
    ----------
    argv:
        Optional sequence of command-line arguments.

    Returns
    -------
    int
        Exit code returned by the script.
    """

    return _run("get_activity_data", argv)


__all__ = [
    "ActivityCommandOptions",
    "MIN_ACTIVITY_TIMEOUT",
    "run_activity_pipeline",
    "main",
]

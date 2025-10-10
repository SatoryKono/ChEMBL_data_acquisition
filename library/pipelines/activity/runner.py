from __future__ import annotations

import argparse
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter

from ... import io
from ...common.log import logger
from ...config import Config
from ...pipelines.assay.chembl_assay import MAX_ACTIVITY_CHUNK_SIZE

MIN_ACTIVITY_TIMEOUT = 60.0

_ACTIVITY_RUNNER: Callable[[Config, argparse.Namespace], int] | None = None
_ACTIVITY_COMPLETION: Callable[..., None] | None = None


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

    if value is None or value is argparse.SUPPRESS:
        return None
    if isinstance(value, Path):
        return value
    return Path(value)


def _normalise_path(value: Path | str) -> Path:
    """Convert a mandatory path-like value to :class:`Path`."""

    result = _normalise_optional_path(value)
    if result is None:  # pragma: no cover - defensive programming
        msg = "a concrete path value is required"
        raise ValueError(msg)
    return result


def register_activity_pipeline_hooks(
    *,
    runner: Callable[[Config, argparse.Namespace], int],
    emit_completion_message: Callable[..., None],
) -> None:
    """Register default dependencies for the activity pipeline runner."""

    global _ACTIVITY_RUNNER, _ACTIVITY_COMPLETION
    _ACTIVITY_RUNNER = runner
    _ACTIVITY_COMPLETION = emit_completion_message


def resolve_activity_pipeline_hooks() -> tuple[
    Callable[[Config, argparse.Namespace], int],
    Callable[..., None],
]:
    """Return the registered activity pipeline dependencies."""

    global _ACTIVITY_RUNNER, _ACTIVITY_COMPLETION

    if _ACTIVITY_RUNNER is None or _ACTIVITY_COMPLETION is None:
        _register_default_activity_pipeline_hooks()

    if _ACTIVITY_RUNNER is None or _ACTIVITY_COMPLETION is None:
        try:  # pragma: no cover - defensive fallback for script import side-effects
            from scripts import get_activity_data as activity_cli  # type: ignore
        except Exception:  # pragma: no cover - the runner may be registered elsewhere
            activity_cli = None
        else:
            register_activity_pipeline_hooks(
                runner=activity_cli.run_chembl,
                emit_completion_message=activity_cli._emit_completion_message,
            )

    if _ACTIVITY_RUNNER is None or _ACTIVITY_COMPLETION is None:
        msg = "activity pipeline runner dependencies could not be resolved"
        raise RuntimeError(msg)

    return _ACTIVITY_RUNNER, _ACTIVITY_COMPLETION


def _register_default_activity_pipeline_hooks() -> None:
    """Register the default activity pipeline handlers from the CLI command."""

    global _ACTIVITY_RUNNER, _ACTIVITY_COMPLETION

    if _ACTIVITY_RUNNER is not None and _ACTIVITY_COMPLETION is not None:
        return

    try:
        from ...cli.commands import get_activity_data as activity_cli
    except Exception:  # pragma: no cover - defensive fallback to scripts
        return

    runner = getattr(activity_cli, "run_chembl", None)
    emit_completion = getattr(activity_cli, "_emit_completion_message", None)
    if runner is None or emit_completion is None:
        return

    register_activity_pipeline_hooks(
        runner=runner,
        emit_completion_message=emit_completion,
    )


_register_default_activity_pipeline_hooks()


def run_activity_pipeline(
    cfg: Config,
    options: ActivityCommandOptions,
    *,
    runner: Callable[[Config, argparse.Namespace], int] | None = None,
    emit_completion_message: Callable[..., None] | None = None,
) -> int:
    """Execute the activity pipeline with orchestrator semantics."""

    if runner is None or emit_completion_message is None:
        default_runner, default_completion = resolve_activity_pipeline_hooks()
        runner = runner or default_runner
        emit_completion_message = emit_completion_message or default_completion

    start_time = perf_counter()
    effective_cfg = cfg.model_copy(deep=True)

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
        else getattr(getattr(effective_cfg, "activity", object()), "timeout", None)
    )
    args.batch_size = (
        options.batch_size
        if options.batch_size is not None
        else getattr(getattr(effective_cfg, "activity", object()), "batch_size", None)
    )
    args.workers = (
        options.workers
        if options.workers is not None
        else getattr(getattr(effective_cfg, "activity", object()), "workers", None)
    )
    if options.invocation is not None:
        args.invocation = tuple(str(part) for part in options.invocation)

    final_output = _normalise_optional_path(options.final_output)
    legacy_output = _normalise_optional_path(options.output_csv)
    if final_output is None:
        if legacy_output is not None:
            output_path = legacy_output
        else:
            output_path = Path(io.default_output_path(args.input_csv, effective_cfg.io))
        args.final_out = output_path
        args.output_csv = output_path
    else:
        output_path = final_output
        args.final_out = output_path
        args.output_csv = output_path

    batch_size = getattr(effective_cfg.activity, "batch_size", None)
    if batch_size is not None and batch_size > MAX_ACTIVITY_CHUNK_SIZE:
        logger.warning(
            "activity_batch_size_clamped",
            configured=batch_size,
            limit=MAX_ACTIVITY_CHUNK_SIZE,
        )
        logger.warning(
            "Configured batch size %s exceeds the hard cap of %s; reducing to %s.",
            batch_size,
            MAX_ACTIVITY_CHUNK_SIZE,
            MAX_ACTIVITY_CHUNK_SIZE,
        )
        effective_cfg.activity.batch_size = MAX_ACTIVITY_CHUNK_SIZE
        args.batch_size = MAX_ACTIVITY_CHUNK_SIZE

    timeout = getattr(effective_cfg.activity, "timeout", None)
    if timeout is not None and timeout < MIN_ACTIVITY_TIMEOUT:
        logger.warning(
            "activity_timeout_clamped",
            configured=timeout,
            minimum=MIN_ACTIVITY_TIMEOUT,
        )
        logger.warning(
            "Configured timeout %s is below the minimum of %s; increasing to %s.",
            timeout,
            MIN_ACTIVITY_TIMEOUT,
            MIN_ACTIVITY_TIMEOUT,
        )
        minimum_timeout = float(MIN_ACTIVITY_TIMEOUT)
        effective_cfg.activity.timeout = minimum_timeout
        try:
            args.timeout = float(minimum_timeout)
        except (TypeError, ValueError):  # pragma: no cover - defensive
            args.timeout = minimum_timeout

    retry_attempts = getattr(effective_cfg.retry, "max_attempts", None)
    if retry_attempts is not None and retry_attempts <= 1:
        logger.warning(
            "activity_retry_disabled",
            configured=retry_attempts,
        )
        logger.warning(
            "Configured system.retry.max_attempts=%s disables urllib3 retry handling; "
            "increase to at least 2 to tolerate transient network issues.",
            retry_attempts,
        )

    api_retries = getattr(effective_cfg.api, "retries", None)
    if api_retries is not None and api_retries <= 0:
        logger.warning(
            "activity_api_retry_disabled",
            configured=api_retries,
        )
        logger.warning(
            "Configured chembl.api.retries=%s disables client-level request retries; "
            "increase to 1 or more for resilience.",
            api_retries,
        )

    if args.skip_existing and output_path.exists() and not args.force:
        logger.info("pipeline_skip_existing", output=str(output_path))
        logger.info(
            "pipeline_skip_existing_reason",
            output=str(output_path),
            message=(
                "Skipping execution because the target output already exists and "
                "--force was not provided."
            ),
        )
        emit_completion_message(
            output_path=output_path,
            processed_rows=None,
            duration_s=perf_counter() - start_time,
            mode="skip_existing",
        )
        return 0

    return runner(effective_cfg, args)


__all__ = [
    "ActivityCommandOptions",
    "MIN_ACTIVITY_TIMEOUT",
    "register_activity_pipeline_hooks",
    "resolve_activity_pipeline_hooks",
    "run_activity_pipeline",
]

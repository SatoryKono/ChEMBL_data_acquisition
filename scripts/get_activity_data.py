"""Compatibility wrapper exposing the activity pipeline CLI entry point."""

from __future__ import annotations

# ruff: noqa: E402  # bootstrap alters import order for script compatibility
import argparse
import json
import sys
from dataclasses import dataclass
from functools import partial, wraps
from importlib import import_module
from itertools import islice
from pathlib import Path
from threading import Lock
from time import perf_counter, sleep
from typing import Any, Callable, Iterable, Iterator, Mapping, Sequence, cast
from urllib.parse import urlsplit

if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli
else:  # pragma: no cover - executed when imported as a module
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli

from library.cli.entrypoints import activity as _activity
from library.pipelines.activity import runner as _activity_runner
from library.pipelines.activity.runner import register_activity_pipeline_hooks
from library.pipelines.assay.chembl_assay import MAX_ACTIVITY_CHUNK_SIZE

ActivityCommandOptions = _activity.ActivityCommandOptions
ActivityPipelineCLI = _activity.ActivityPipelineCLI
PreparedActivityContext = _activity.PreparedActivityContext

DEFAULT_INPUT_NAME = _activity.DEFAULT_INPUT_NAME
DEFAULT_OUTPUT_STEM = _activity.DEFAULT_OUTPUT_STEM
PROGRAM_NAME = _activity.PROGRAM_NAME
MIN_ACTIVITY_TIMEOUT = _activity.MIN_ACTIVITY_TIMEOUT

cl = _activity.cl
configure_logger = _activity.configure_logger
file_sha256 = _activity.file_sha256
io = _activity.io
logger = _activity.logger
run_cli_command = _activity.run_cli_command
sleep = _activity.sleep
write_csv_chunks_deterministic = _activity.write_csv_chunks_deterministic
write_meta_yaml = _activity.write_meta_yaml

_args_invocation = _activity._args_invocation


def _sync_runtime_dependencies() -> None:
    module_globals = globals()
    logger_obj = module_globals["logger"]
    _activity.logger = logger_obj
    _activity_runner.logger = logger_obj
    _activity.sleep = module_globals["sleep"]
    _activity.write_csv_chunks_deterministic = module_globals[
        "write_csv_chunks_deterministic"
    ]
    _activity.run_chembl = module_globals["run_chembl"]
    _activity._emit_completion_message = module_globals["_emit_completion_message"]


def _with_runtime_sync(func):
    from functools import wraps
    @wraps(func)
    def wrapper(*args, **kwargs):
        _sync_runtime_dependencies()
        return func(*args, **kwargs)
    return wrapper
    ]
    _activity.run_chembl = module_globals["run_chembl"]
    _activity._emit_completion_message = module_globals["_emit_completion_message"]


def _with_runtime_sync(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        _sync_runtime_dependencies()
        return func(*args, **kwargs)

    return wrapper

    last_error_extra: dict[str, object] | None = None
    last_error_context: dict[str, object] = {}

    pref_name_cache: dict[str, str | None] = {}
    pref_name_cache_lock = Lock()

    with ETLContext(cfg) as context:
        client = context.chembl_client

        retry_cfg = cfg.retry
        chunk_failures = ChunkFailureTracker()
        http_diagnostics = _gather_http_diagnostics(cfg, client)
        if http_diagnostics:
            logger.info("activity_http_config", **http_diagnostics)

        def fetch_chunk(chunk_ids: Sequence[str]) -> pd.DataFrame:
            nonlocal last_error_extra, last_error_context

            timeout_error_types = (
                requests.Timeout,
                requests.ReadTimeout,
                requests.ConnectTimeout,
                requests.exceptions.RetryError,
            )

            def _is_timeout_error(exc: Exception) -> bool:
                if isinstance(exc, timeout_error_types):
                    return True
                message = str(exc).strip().lower()
                if not message:
                    return False
                return "timed out" in message or "timeout" in message

            def _fetch(ids: Sequence[str], *, depth: int = 0) -> pd.DataFrame:
                nonlocal last_error_extra, last_error_context

                attempts = max(1, retry_cfg.max_attempts)
                id_list = [str(identifier) for identifier in ids]

                for attempt in range(1, attempts + 1):
                    try:
                        result = cl.get_activities(
                            id_list,
                            cfg=cfg.api,
                            client=client,
                            chunk_size=cfg.activity.batch_size,
                            timeout=cfg.activity.timeout,
                            **extra_kwargs,
                        )
                    except (requests.RequestException, ValueError) as exc:
                        error_message = str(exc)
                        network_hint, network_host = _describe_network_failure(cfg, exc)
                        context = {
                            "chunk_ids": list(id_list),
                            "chunk_size": len(id_list),
                            "attempt": attempt,
                            "max_attempts": attempts,
                            "batch_size": cfg.activity.batch_size,
                            "timeout": cfg.activity.timeout,
                        }
                        if network_host:
                            context["api_host"] = network_host
                        for key in (
                            "adapter_total",
                            "api_retries",
                            "retry_max_attempts",
                            "activity_timeout",
                            "api_timeout_read",
                        ):
                            if key in http_diagnostics:
                                context[key] = http_diagnostics[key]
                        log_context = {
                            key: value for key, value in context.items() if key != "chunk_ids"
                        }
                        last_error_extra = {
                            "msg": error_message,
                            "chunk_ids": context["chunk_ids"],
                        }
                        if network_hint:
                            last_error_extra["hint"] = network_hint
                        last_error_context = dict(log_context)
                        if network_host:
                            last_error_context.setdefault("api_host", network_host)
                        if network_hint:
                            last_error_context["network_hint"] = network_hint
                        if attempt >= attempts:
                            if network_hint:
                                logger.error(
                                    "activity_fetch_network_error",
                                    extra=last_error_extra,
                                    hint=network_hint,
                                    **log_context,
                                )
                            if len(id_list) > 1 and _is_timeout_error(exc):
                                split_context = dict(log_context)
                                split_context["depth"] = depth
                                logger.warning(
                                    "activity_fetch_split",
                                    extra=last_error_extra,
                                    **split_context,
                                )
                                midpoint = max(1, len(id_list) // 2)
                                left_ids = tuple(id_list[:midpoint])
                                right_ids = tuple(id_list[midpoint:])
                                frames: list[pd.DataFrame] = []
                                if left_ids:
                                    frames.append(_fetch(left_ids, depth=depth + 1))
                                if right_ids:
                                    frames.append(_fetch(right_ids, depth=depth + 1))
                                if frames:
                                    combined = pd.concat(
                                        frames, ignore_index=True, sort=False
                                    )
                                else:
                                    combined = pd.DataFrame(columns=ACTIVITY_COLUMNS)
                                last_error_extra = None
                                last_error_context = {}
                                return combined
                            logger.error(
                                "activity_fetch_failed",
                                extra=last_error_extra,
                                error=error_message,
                                **log_context,
                            )
                            chunk_failures.add_failure(id_list, error_message)
                            raise PipelineError("chunk_fetch_failed")  # noqa: B904
                        delay = compute_backoff_delay(attempt, retry_cfg)
                        logger.warning(
                            "activity_fetch_retry",
                            extra=last_error_extra,
                            delay=delay,
                            **log_context,
                        )
                        if delay > 0:
                            sleep(delay)
                    else:
                        last_error_extra = None
                        last_error_context = {}
                        return _ensure_molecule_pref_name(
                            result,
                            cfg=cfg,
                            client=client,
                            cache=pref_name_cache,
                            cache_lock=pref_name_cache_lock,
                            chunk_failures=chunk_failures,
                        )
                return pd.DataFrame(columns=ACTIVITY_COLUMNS)

            return _fetch(chunk_ids, depth=0)

        worker_count = getattr(cfg.activity, "workers", 1) or 1
        fetch_config = ChunkedFetchConfig(
            ids=limited_ids,
            chunk_size=cfg.activity.batch_size,
            workers=max(1, worker_count),
        )

        writer_config = CsvWriterConfig(
            writer=writer,
            kwargs={},
            ensure_destination=True,

        )

        pipeline_stats: dict[str, object] | None = None

        def _capture_stats(stats: Mapping[str, object]) -> None:
            nonlocal pipeline_stats, streamed_summary
            pipeline_stats = dict(stats)
            snapshot = streaming_stats.snapshot()
            streamed_summary = snapshot
            pipeline_stats.setdefault("rows_streamed", snapshot["rows"])
            pipeline_stats.setdefault("cells_streamed", snapshot["cells"])
            pipeline_stats.setdefault("null_cells", snapshot["nulls"])
            pipeline_stats.setdefault("null_fraction", snapshot["null_fraction"])

        definition_kwargs: dict[str, object] = {
            "schema": ActivitiesSchema,
            "schema_name": "ActivitiesSchema",
            "validators": validators,
            "command": " ".join(_args_invocation(args)),
            "config_snapshot": _serialize_paths(cfg.to_dict()),
            "inputs": {"input_csv": str(args.input_csv)},
            "key_columns": ["activity_id"],
            "table_quality": table_quality,
            "stats_extra": chunk_failures.stats,
            "stats_callback": _capture_stats,
            "dictionary_resources": (
                "dictionary_root",
                "target_types",
            ),
        }

        result = activity_run.run_activity_pipeline(
            fetch_config=fetch_config,
            metadata_hooks=metadata_hooks,
            fetch_chunk=fetch_chunk,
            writer_config=writer_config,
            definition_kwargs=definition_kwargs,
            cfg=cfg,
            logger=logger,
            output_path=output_path,
            failure_path=failure_path,
            fetch_failure_path=fetch_failure_path,
            chunk_failures=chunk_failures,
        )
        exit_code = result.exit_code

    if exit_code == 0:
        logger.info(
            f"Merged data checkpoint: wrote merged activity records to '{output_path}'."
        )
        try:
            extended_output_path = process_activity_extended(
                input_path=output_path,
                search_dir=output_path.parent,
                dictionary_dir=cfg.resources.dictionary_dir,
                targets_csv=cfg.resources.targets_type_csv,
            )
        except Exception:
            logger.exception(
                "Activity extended post-processing failed while enriching merged activity data."
            )
            raise

    processed_ids = prepared_context.processed_ids
    processed_count = 0
    try:
        processed_count = int(processed_ids or 0)
    except (TypeError, ValueError):
        logger.info("processed_count_conversion_failed", value=processed_ids)
        processed_count = 0

    summary_snapshot: dict[str, object] | None = None

    if limit is not None:
        logger.info(
            "process_limit",
            processed=processed_ids,
            limit=limit,
        )
        summary_snapshot = streamed_summary or streaming_stats.snapshot()
    else:
        summary_snapshot = streamed_summary or streaming_stats.snapshot()
        logger.info("processed_count", count=processed_count)

    if pipeline_stats is not None:
        try:
            rows_total = int(pipeline_stats.get("rows_total", processed_count))
        except (TypeError, ValueError):
            rows_total = processed_count
        try:
            rows_kept = int(pipeline_stats.get("rows_kept", 0))
        except (TypeError, ValueError):
            rows_kept = 0
        rows_dropped = int(pipeline_stats.get("rows_dropped", 0))
        logger.info(
            "records_dropped",
            rows_total=rows_total,
            rows_kept=rows_kept,
            rows_dropped=rows_dropped,
        )

    if exit_code == 0:
        completion_rows = processed_ids
        summary_rows = summary_snapshot.get("rows") if summary_snapshot else None
        if isinstance(summary_rows, int | float):
            try:
                completion_rows = int(summary_rows)
            except Exception:
                # Fallback in case of any unexpected type or value
                completion_rows = processed_ids
        elif pipeline_stats is not None:
            try:
                # Try to get from pipeline_stats if available
                rows_kept = pipeline_stats.get("rows_kept")
                if rows_kept is not None and isinstance(rows_kept, int | float):
                    completion_rows = int(rows_kept)
                else:
                    rows_total = pipeline_stats.get("rows_total", processed_ids)
                    if isinstance(rows_total, int | float):
                        completion_rows = int(rows_total)
                    else:
                        completion_rows = processed_ids
            except Exception:
                # Fallback in case of any unexpected type or value
                completion_rows = processed_ids
        else:
            completion_rows = processed_ids

        report_extras: dict[str, object] = {"rows": completion_rows, "processed": processed_ids}
        if pipeline_stats is not None:
            report_extras.update(pipeline_stats)
        if summary_snapshot:
            report_extras["summary_snapshot"] = summary_snapshot
        if extended_output_path is not None:
            report_extras["extended_output"] = str(extended_output_path)

        postprocess_metrics, report_path = _generate_activity_postprocess_metrics(
            cfg,
            output_path,
            logger=logger,
            extras=report_extras,
        )

        pipeline_version_value = (
            getattr(postprocess_metrics, "pipeline_version", None)
            if postprocess_metrics and postprocess_metrics.pipeline_version is not None
            else get_pipeline_version()
        )

        pipeline_done_payload: dict[str, object] = {
            "output": str(output_path),
            "rows": completion_rows,
            "pipeline_version": pipeline_version_value,
        }
        if summary_snapshot is not None:
            null_fraction = summary_snapshot.get("null_fraction")
            if null_fraction is not None:
                pipeline_done_payload["null_fraction"] = null_fraction
        if extended_output_path is not None:
            pipeline_done_payload["extended_output"] = str(extended_output_path)
        if postprocess_metrics is not None:
            summary = postprocess_metrics.summary()
            if summary.get("rows") is not None:
                pipeline_done_payload["postprocess_rows"] = summary["rows"]
            if summary.get("columns") is not None:
                pipeline_done_payload["postprocess_columns"] = summary["columns"]
            if summary.get("duration_s") is not None:
                pipeline_done_payload["postprocess_duration_s"] = summary["duration_s"]
            if summary.get("steps") is not None:
                pipeline_done_payload["postprocess_steps"] = summary["steps"]
            validation = getattr(postprocess_metrics, "validation", None)
            if validation is not None:
                pipeline_done_payload["postprocess_schema"] = getattr(validation, "schema", None)
            if report_path is not None:
                pipeline_done_payload["postprocess_report"] = str(report_path)

        logger.info("activity_pipeline_done", extra=pipeline_done_payload)
        if extended_output_path is not None:
            logger.info(
                "activity_export_checkpoint",
                output=str(output_path),
                extended_output=str(extended_output_path),
            )
        else:
            logger.info(
                "activity_export_checkpoint",
                output=str(output_path),
                extended_output=None,
            )
        _emit_completion_message(
            output_path=output_path,
            processed_rows=completion_rows,
            duration_s=perf_counter() - start_time,
            mode="run",
            streamed_metrics=summary_snapshot,
        )
    else:
        extra_payload = last_error_extra
        context_payload = dict(last_error_context) if last_error_context else {}
        error_message = None
        chunk_ids: Sequence[str] | None = None
        if extra_payload:
            error_message = extra_payload.get("msg")
            chunk_ids = extra_payload.get("chunk_ids")  # type: ignore[assignment]
        details = []
        if error_message:
            details.append(f"last error: {error_message}")
        if chunk_ids:
            details.append(f"chunk_ids={list(chunk_ids)}")
        if context_payload:
            attempt_info = ", ".join(f"{key}={value}" for key, value in context_payload.items())
            details.append(attempt_info)
        detail_text = "; ".join(details)
        logger.error(
            "activity_pipeline_failed",
            output=str(output_path),
            processed=processed_ids,
            exit_code=exit_code,
            details=detail_text or None,
        )
        logger.error(
            "activity_pipeline_failed_detail",
            output=str(output_path),
            processed=processed_ids,
            exit_code=exit_code,
            details=detail_text or None,
        )

    return exit_code

_emit_completion_message = _with_runtime_sync(_activity._emit_completion_message)
_ensure_extended_activity_columns = _with_runtime_sync(
    _activity._ensure_extended_activity_columns
)
_ensure_molecule_pref_name = _with_runtime_sync(_activity._ensure_molecule_pref_name)
_ensure_src_assay_id = _with_runtime_sync(_activity._ensure_src_assay_id)
_filter_activity_output_columns = _with_runtime_sync(
    _activity._filter_activity_output_columns
)
_load_assay_src_lookup = _with_runtime_sync(_activity._load_assay_src_lookup)
main = _with_runtime_sync(_activity.main)
prepare_activity_context = _with_runtime_sync(_activity.prepare_activity_context)
run = _with_runtime_sync(_activity.run)
run_chembl = _with_runtime_sync(_activity.run_chembl)

register_activity_pipeline_hooks(
    runner=run_chembl,
    emit_completion_message=_emit_completion_message,
)

__all__ = (
    "ActivityCommandOptions",
    "ActivityPipelineCLI",
    "DEFAULT_INPUT_NAME",
    "DEFAULT_OUTPUT_STEM",
    "MAX_ACTIVITY_CHUNK_SIZE",
    "MIN_ACTIVITY_TIMEOUT",
    "PROGRAM_NAME",
    "PreparedActivityContext",
    "cl",
    "configure_logger",
    "file_sha256",
    "io",
    "logger",
    "main",
    "prepare_activity_context",
    "run",
    "run_chembl",
    "run_cli_command",
    "sleep",
    "write_csv_chunks_deterministic",
    "write_meta_yaml",
    "_args_invocation",
    "_emit_completion_message",
    "_ensure_extended_activity_columns",
    "_ensure_molecule_pref_name",
    "_ensure_src_assay_id",
    "_filter_activity_output_columns",
    "_load_assay_src_lookup",
)


def __getattr__(name: str) -> object:
    """Proxy attribute access to :mod:`library.cli.entrypoints.activity`."""

    try:
        return getattr(_activity, name)
    except AttributeError as exc:  # pragma: no cover - passthrough for missing attrs
        raise AttributeError(name) from exc


from library.cli.entrypoints import activity as _activity

sys.modules[__name__] = _activity
sys.modules.setdefault("scripts.get_activity_data", _activity)

if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(_activity.main())

"""Helpers shared across the command line interfaces.

This module exposes utilities for constructing ``argparse`` parsers and
coordinating the execution of data extraction pipelines.  The
``run_pipeline`` helper consolidates the boilerplate shared by multiple
scripts: fetching raw data in chunks, applying metadata hooks, validating
against Pandera schemas and writing deterministic CSV outputs alongside
metadata files.
"""

from __future__ import annotations

import argparse
import shlex
import sys
import traceback
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from contextvars import ContextVar
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
from pandera.errors import SchemaErrors

from ..cli import (
    Logger,
    LoggerConfig,
    add_common_arguments,
    apply_config_overrides,
    configure_logger,
    path_argument,
    positive_int,
    prepare_io_paths,
)
from ..cli_utils import ensure_run_id
from ..common.log import logger as default_logger
from ..common.metadata import record_quality_failure
from ..common.sidecar import SidecarErrors
from ..config import (
    DEFAULT_CONFIG_PATH,
    Config,
    ConfigError,
    ensure_dirs,
    print_config,
)
from ..maintenance import ensure_legacy_cleanup

if TYPE_CHECKING:
    from ..postprocess import PostprocessingPipelineResult
from ..reporting.run_manifest import finalise_csv_output
from .pipeline_definition import (
    Fetcher,
    PipelineDefinition,
    normalise_definition,
)


@dataclass(frozen=True)
class _PostprocessHandlers:
    runner: Callable[..., tuple[pd.DataFrame, object]]
    validator: Callable[..., pd.DataFrame]
    schema: object


_POSTPROCESS_TABLE_ALIASES = {
    "activity": "activities",
    "activities": "activities",
    "assay": "assays",
    "assays": "assays",
    "document": "documents",
    "documents": "documents",
    "target": "targets",
    "targets": "targets",
    "testitem": "testitem",
    "testitems": "testitem",
}


@dataclass(frozen=True)
class _PostprocessRuntime:
    pipeline_config_factory: Callable[..., object]
    run_pipeline: Callable[..., object]
    load_pipeline_config: Callable[[str, Path | str | None], object]
    get_csv_runtime_config: Callable[[object], object]
    handlers: Mapping[str, _PostprocessHandlers]


@dataclass(frozen=True)
class _PostprocessOptions:
    enabled: bool
    config_override: Path | str | None


_POSTPROCESS_OPTIONS: ContextVar[_PostprocessOptions | None] = ContextVar(
    "_POSTPROCESS_OPTIONS", default=None
)


def _get_postprocess_options() -> _PostprocessOptions:
    options = _POSTPROCESS_OPTIONS.get(None)
    if options is None:
        return _PostprocessOptions(False, None)
    return options


@lru_cache(maxsize=1)
def _load_postprocess_runtime() -> _PostprocessRuntime:
    from ..postprocess.common import (
        PostprocessingPipelineConfig,
        run_postprocessing_pipeline,
    )
    from ..postprocess.common import (
        get_csv_runtime_config as get_postprocess_csv_config,
    )
    from ..postprocess.common import (
        get_pipeline_config as load_postprocess_pipeline_config,
    )
    from ..postprocessing.activities import (
        ACTIVITY_SCHEMA,
        validate_activities,
    )
    from ..postprocessing.activities import (
        run_activity_pipeline as run_activity_postprocess,
    )
    from ..postprocessing.assays import (
        ASSAY_SCHEMA,
        validate_assays,
    )
    from ..postprocessing.assays import (
        run_assay_pipeline as run_assay_postprocess,
    )
    from ..postprocessing.documents import (
        DOCUMENT_SCHEMA,
        validate_documents,
    )
    from ..postprocessing.documents import (
        run_document_pipeline as run_document_postprocess,
    )
    from ..postprocessing.targets import (
        TARGET_SCHEMA,
        validate_targets,
    )
    from ..postprocessing.targets import (
        run_target_pipeline as run_target_postprocess,
    )
    from ..postprocessing.testitem import (
        TESTITEM_SCHEMA,
        validate_testitems,
    )
    from ..postprocessing.testitem import (
        run_testitem_pipeline as run_testitem_postprocess,
    )

    handlers: dict[str, _PostprocessHandlers] = {
        "activities": _PostprocessHandlers(
            runner=run_activity_postprocess,
            validator=validate_activities,
            schema=ACTIVITY_SCHEMA,
        ),
        "assays": _PostprocessHandlers(
            runner=run_assay_postprocess,
            validator=validate_assays,
            schema=ASSAY_SCHEMA,
        ),
        "documents": _PostprocessHandlers(
            runner=run_document_postprocess,
            validator=validate_documents,
            schema=DOCUMENT_SCHEMA,
        ),
        "targets": _PostprocessHandlers(
            runner=run_target_postprocess,
            validator=validate_targets,
            schema=TARGET_SCHEMA,
        ),
        "testitem": _PostprocessHandlers(
            runner=run_testitem_postprocess,
            validator=validate_testitems,
            schema=TESTITEM_SCHEMA,
        ),
        "testitems": _PostprocessHandlers(
            runner=run_testitem_postprocess,
            validator=validate_testitems,
            schema=TESTITEM_SCHEMA,
        ),
    }

    return _PostprocessRuntime(
        pipeline_config_factory=PostprocessingPipelineConfig,
        run_pipeline=run_postprocessing_pipeline,
        load_pipeline_config=load_postprocess_pipeline_config,
        get_csv_runtime_config=get_postprocess_csv_config,
        handlers=handlers,
    )


def _extract_output_token(path: Path) -> str | None:
    name = path.name
    if name.startswith("output.") and name.endswith(".csv"):
        return name[len("output.") : -len(".csv")]
    return None


def _resolve_postprocess_table(token: str) -> str | None:
    base = token.split("_", 1)[0].lower().strip()
    return _POSTPROCESS_TABLE_ALIASES.get(base)


def _maybe_run_postprocessing(
    csv_path: Path,
    logger: Logger,
    *,
    postprocess_enabled: bool,
    config_override: Path | str | None,
) -> PostprocessingPipelineResult | None:
    token = _extract_output_token(csv_path)
    if not token:
        return None

    table = _resolve_postprocess_table(token)
    if table is None:
        return None

    if not postprocess_enabled:
        logger.info("[INFO] Postprocessing skipped (flag --postprocess not set)")
        return None

    try:
        runtime = _load_postprocess_runtime()
        handlers = runtime.handlers.get(table)
        if handlers is None:
            return None
        pipeline_cfg = runtime.load_pipeline_config(table, config_override)
        csv_cfg = runtime.get_csv_runtime_config(pipeline_cfg)
        pipeline_config = runtime.pipeline_config_factory(
            pipeline_config=pipeline_cfg,
            csv_runtime_config=csv_cfg,
            runner=handlers.runner,
            validator=handlers.validator,
            schema=handlers.schema,
            logger=logger,
        )
        destination = csv_path.with_name(f"output_postprocessed.{token}.csv")
        result = runtime.run_pipeline(
            table,
            csv_path,
            destination,
            pipeline_config,
        )
    except Exception:
        logger.exception(
            "postprocess_failed",
            table=table,
            input=str(csv_path),
        )
        raise
    else:
        logger.info(
            "postprocess_done",
            table=table,
            input=str(csv_path),
            output=str(result.output_path),
        )
        return result


def _callable_name(func: Callable[..., object]) -> str:
    """Return a human readable name for ``func``."""

    return getattr(func, "__qualname__", getattr(func, "__name__", repr(func)))


class PipelineError(RuntimeError):
    """Raised when a pipeline step encounters a fatal error."""


def resolve_invocation(
    program: str | None, argv: Sequence[object] | None
) -> tuple[str, ...]:
    """Return a normalised tuple describing the CLI invocation.

    Parameters
    ----------
    program:
        Optional executable name. When provided it is prefixed to the
        resulting tuple.
    argv:
        Sequence of arguments. ``None`` falls back to ``sys.argv[1:]`` to
        mirror standard ``argparse`` semantics.

    Returns
    -------
    tuple[str, ...]
        Normalised invocation containing the program followed by positional
        arguments represented as strings.
    """

    parts: list[str] = []
    if program:
        parts.append(str(program))
    if argv is None:
        argv = sys.argv[1:]
    parts.extend(str(arg) for arg in argv)
    return tuple(parts)


def _normalise_path(value: object) -> str | None:
    if isinstance(value, Path):
        candidate = value
    elif isinstance(value, str):
        candidate = Path(value)
    elif hasattr(value, "__fspath__"):
        candidate = Path(value)
    else:
        return None
    try:
        return str(candidate.expanduser().resolve())
    except OSError:
        return str(candidate.expanduser())


def _canonical_run_descriptor(
    args: argparse.Namespace, parser: argparse.ArgumentParser
) -> str:
    invocation = getattr(args, "invocation", resolve_invocation(parser.prog, None))
    parts = [str(part) for part in invocation]
    path_fields = (
        "base_path",
        "input_dir",
        "output_dir",
        "cache_dir",
        "input_csv",
        "output_csv",
        "final_out",
        "raw_out",
    )
    entries: list[str] = []
    for field in path_fields:
        value = getattr(args, field, None)
        if value in (None, argparse.SUPPRESS):
            continue
        normalised = _normalise_path(value)
        if normalised is not None:
            entries.append(f"{field}={normalised}")
    if entries:
        parts.extend(sorted(entries))
    return "\n".join(parts)


def run_cli_command(
    *,
    args: argparse.Namespace,
    parser: argparse.ArgumentParser,
    base_parser: argparse.ArgumentParser | None = None,
    log_cfg: LoggerConfig,
    mapping: Mapping[str, str],
    run: Callable[[Config, argparse.Namespace], int],
    logger: Logger | None = None,
    postprocess_enabled: bool | None = None,
) -> int:
    """Execute CLI boilerplate shared by data acquisition commands."""

    debug_flag = bool(getattr(args, "debug", False))
    keep_flag = bool(getattr(args, "keep_intermediate", False))
    legacy_flag = bool(getattr(args, "emit_legacy_artifacts", False))
    diagnostics_enabled = legacy_flag or debug_flag or keep_flag
    if diagnostics_enabled != legacy_flag and hasattr(args, "emit_legacy_artifacts"):
        args.emit_legacy_artifacts = diagnostics_enabled

    level_candidate = getattr(args, "log_level", log_cfg.level)
    level = str(level_candidate).upper()
    if getattr(args, "verbose", False):
        level = "DEBUG"
    log_cfg.level = level

    def _configure_logging() -> Logger:
        configured_logger = configure_logger(log_cfg)
        return logger or configured_logger

    try:
        config_arg = getattr(args, "config", None)
        if isinstance(config_arg, str | Path):
            config_path: Path | str = config_arg
        else:
            default_config = parser.get_default("config")
            if not isinstance(default_config, str | Path):
                msg = "configuration path must be provided"
                raise ValueError(msg)
            config_path = default_config
    except ValueError as exc:
        ensure_run_id(args, parser, log_cfg)
        use_logger = _configure_logging()
        use_logger.error(
            "config_error",
            error=str(exc),
            config=str(getattr(args, "config", "")),
            exc_info=exc,
        )
        use_logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1

    try:
        prepare_io_paths(args)
    except (ValueError, FileNotFoundError) as exc:
        ensure_run_id(args, parser, log_cfg)
        use_logger = _configure_logging()
        use_logger.error(
            "config_error",
            error=str(exc),
            config=str(config_path),
            exc_info=exc,
        )
        use_logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1

    ensure_run_id(args, parser, log_cfg)
    use_logger = _configure_logging()

    use_logger.info("pipeline_start", run_id=log_cfg.run_id)

    try:
        cfg: Config = apply_config_overrides(
            args,
            parser,
            config_path,
            mapping=dict(mapping),
            base_parser=base_parser,
        )
        metadata = getattr(args, "_config_metadata", None)
    except (ConfigError, FileNotFoundError, ValueError) as exc:
        use_logger.error(
            "config_error",
            error=str(exc),
            config=str(config_path),
            exc_info=exc,
        )
        use_logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1

    if metadata is not None:
        use_logger.info("config_snapshot", config=getattr(metadata, "snapshot", {}))

    if not diagnostics_enabled:
        system_cfg = getattr(cfg, "system", None)
        doc_quality_cfg = getattr(system_cfg, "doc_quality", None)
        if doc_quality_cfg is not None and hasattr(doc_quality_cfg, "enable"):
            doc_quality_cfg.enable = False

    try:
        if getattr(args, "print_config", False):
            print_config(cfg)
            configure_logger(log_cfg)
            use_logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        cleanup_result = ensure_legacy_cleanup(cfg, logger=use_logger)
        if cleanup_result.performed and not cleanup_result.dry_run:
            log_method = (
                use_logger.warning
                if cleanup_result.removed_count
                else use_logger.info
            )
            log_method(
                "legacy_outputs_retention_notice",
                directory=str(cleanup_result.output_dir),
                removed=cleanup_result.removed_count,
                sentinel=str(cleanup_result.sentinel_path),
                hint=(
                    "Legacy diagnostics are now opt-in via --emit-legacy-artifacts/"
                    "--keep-intermediate/--debug"
                ),
            )
        use_logger = configure_logger(log_cfg)
    except (ValueError, TypeError) as exc:
        use_logger.error(
            "config_error",
            error=str(exc),
            config=str(config_path),
            exc_info=exc,
        )
        use_logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        use_logger.error(
            "directory_setup_failed",
            error=str(exc),
            exc_info=exc,
        )
        use_logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1

    if postprocess_enabled is None:
        postprocess_enabled = bool(getattr(args, "postprocess", False))

    postprocess_config: Path | str | None
    config_attr = getattr(args, "config", None)
    if isinstance(config_attr, str | Path):
        postprocess_config = config_attr
    else:
        postprocess_config = config_path

    token = _POSTPROCESS_OPTIONS.set(
        _PostprocessOptions(
            enabled=postprocess_enabled,
            config_override=postprocess_config,
        )
    )
    try:
        exit_code = run(cfg, args)
    finally:
        _POSTPROCESS_OPTIONS.reset(token)

    if exit_code == 0:
        use_logger.info("pipeline_done", run_id=log_cfg.run_id)
    else:
        use_logger.info("pipeline_fail", run_id=log_cfg.run_id)
    return exit_code


def _as_iterable(
    source: pd.DataFrame | Iterable[pd.DataFrame],
) -> Iterator[pd.DataFrame]:
    """Return ``source`` as an iterator of data frames."""

    if isinstance(source, pd.DataFrame):
        yield source
    else:
        yield from source


def run_pipeline(
    *,
    definition: PipelineDefinition | None = None,
    fetcher: Fetcher,
    output_path: Path,
    failure_path: Path,
    cfg: Config | None = None,
    logger: Logger | None = None,
    postprocess_enabled: bool | None = None,
    postprocess_config: Path | str | None = None,
    **legacy_kwargs: object,
) -> int:
    """Execute a data pipeline and write deterministic CSV output.

    Parameters
    ----------
    definition:
        Pipeline configuration bundle describing schema, validators, metadata
        hooks and writer behaviour. When omitted, the legacy keyword arguments
        from older call-sites are used to construct a temporary definition for
        backwards compatibility.
    fetcher:
        Callable returning an iterable of raw :class:`pandas.DataFrame`
        objects.  Each frame represents a chunk of data retrieved from an
        upstream service.
    output_path:
        Destination ``CSV`` file.
    failure_path:
        Path for persisting validation failure cases.
    cfg:
        Optional application configuration forwarded to sidecar metadata.
    logger:
        Optional logger.  Defaults to :data:`library.common.log.logger` when omitted.
    legacy_kwargs:
        Deprecated keyword arguments matching the historical signature. When
        provided they are converted into a :class:`PipelineDefinition`.

    Returns
    -------
    int
        Zero on success, one when validation fails or ``writer`` raises an
        exception.  ``PipelineError`` exceptions raised by ``fetcher`` are
        converted into a non-zero return code.
    """

    use_logger = logger or default_logger

    definition = normalise_definition(definition, legacy_kwargs)

    context_options = _get_postprocess_options()
    effective_postprocess_enabled = (
        postprocess_enabled
        if postprocess_enabled is not None
        else context_options.enabled
    )
    effective_postprocess_config = (
        postprocess_config
        if postprocess_config is not None
        else context_options.config_override
    )

    schema = definition.schema
    schema_name = definition.schema_name
    metadata_hooks = list(definition.metadata_hooks)
    validators = list(definition.validators)
    writer = definition.writer
    config_snapshot = dict(definition.config_snapshot)
    inputs = dict(definition.inputs)
    key_columns = list(definition.key_columns)
    table_quality = definition.table_quality or (lambda _: None)
    stats_extra = definition.stats_extra
    stats_callback = definition.stats_callback
    strict_mode = bool(definition.strict_mode)
    invocation_tuple = tuple(str(part) for part in (definition.invocation or ()))
    command = definition.command

    # NOTE:
    # ``output_path`` and ``failure_path`` are provided as keyword-only
    # arguments.  In certain execution environments (for instance when the
    # function is reflected and re-bound dynamically) direct access to these
    # parameters can incorrectly raise ``NameError`` during optimisation passes.
    # To shield the pipeline from those edge cases we resolve the values via
    # ``locals()`` before re-binding them to local variables.
    output_path_value = locals().get("output_path")
    failure_path_value = locals().get("failure_path")
    if output_path_value is None or failure_path_value is None:
        raise ValueError("run_pipeline requires 'output_path' and 'failure_path'")
    output_path = Path(output_path_value)
    failure_path = Path(failure_path_value)

    if command is not None:
        command_str = command
    elif invocation_tuple:
        command_str = " ".join(shlex.quote(part) for part in invocation_tuple)
    else:
        raise ValueError("run_pipeline requires either 'command' or 'invocation'")

    # ``definition`` already normalises metadata hooks and validators to tuples
    # so converting to ``list`` above is sufficient for mutation.

    if schema is not None:
        required_cols = {
            name
            for name, column in getattr(schema, "columns", {}).items()
            if column.required
        }
        optional_cols = set(getattr(schema, "columns", {})) - required_cols
    else:
        required_cols = set()
        optional_cols = set()

    errors = SidecarErrors()
    rows_total = 0
    rows_kept = 0
    rows_dropped = 0
    total_failures = 0
    exit_code = 0
    present_columns: set[str] = set()
    missing_required_columns: set[str] = set()
    all_columns: set[str] = set()
    validation_enabled = True
    failed_metadata_hooks: set[str] = set()

    try:
        iterable = _as_iterable(fetcher())
    except PipelineError:
        return 1
    except Exception as exc:  # pragma: no cover - exercised in integration tests
        use_logger.error("fetch_failed", error=str(exc))
        return 1

    class _AbortPipeline(RuntimeError):
        """Internal exception raised to abort processing early."""

        def __init__(self, code: int = 1) -> None:
            super().__init__("pipeline aborted")
            self.code = code

    schema_columns: list[str] | None
    if schema is not None:
        schema_columns = list(getattr(schema, "columns", {}))
    else:
        schema_columns = None

    col_order: list[str] = list(schema_columns or [])
    resolved_keys: list[str] = []

    def _refresh_column_tracking(frame: pd.DataFrame) -> None:
        chunk_columns = set(frame.columns)
        present_columns.update(chunk_columns)
        all_columns.update(chunk_columns)

        if schema_columns is not None:
            head = [column for column in schema_columns if column in all_columns]
            tail = sorted(
                column for column in all_columns if column not in schema_columns
            )
            col_order[:] = head + tail
        else:
            col_order[:] = sorted(all_columns)

        resolved_keys[:] = [column for column in key_columns if column in all_columns]

    def _validated_chunks() -> Iterator[pd.DataFrame]:
        nonlocal rows_total, rows_kept, rows_dropped, total_failures, exit_code
        nonlocal validation_enabled

        chunks_emitted = False
        aborted = False
        chunk_index = 0

        try:
            for chunk in iterable:
                if chunk is None:
                    continue

                chunk_index += 1
                chunk_rows_total = len(chunk)
                rows_total += chunk_rows_total

                for hook in metadata_hooks:
                    hook_name = _callable_name(hook)
                    try:
                        chunk = hook(chunk)
                    except Exception as exc:
                        use_logger.error(
                            "metadata_hook_failed",
                            hook=hook_name,
                            error=str(exc),
                            error_type=exc.__class__.__name__,
                            context="chunk",
                            chunk_index=chunk_index,
                            rows=chunk_rows_total,
                            strict_mode=strict_mode,
                        )
                        failed_metadata_hooks.add(hook_name)
                        if strict_mode:
                            aborted = True
                            raise _AbortPipeline(1) from exc

                _refresh_column_tracking(chunk)

                missing_chunk_required = required_cols - set(chunk.columns)
                if missing_chunk_required:
                    missing_required_columns.update(missing_chunk_required)
                    validation_enabled = False
                    exit_code = 1
                    aborted = True
                    raise _AbortPipeline(exit_code)

                validated_chunk = chunk
                if validation_enabled and validators and not chunk.empty:
                    for validator in validators:
                        try:
                            result = validator(validated_chunk)
                        except SchemaErrors as exc:
                            failure_cases = exc.failure_cases
                            if not failure_cases.empty:
                                for row in failure_cases.to_dict("records"):
                                    errors.add_error(row)
                                total_failures += len(failure_cases)
                                exit_code = 1
                                use_logger.error(
                                    "validation_failed",
                                    failures=len(failure_cases),
                                    path=str(failure_path),
                                )
                            validated_chunk = getattr(
                                exc, "validated_data", validated_chunk
                            )
                        else:
                            validated_chunk = result.data
                            failure_cases = result.failure_cases
                            if not failure_cases.empty:
                                for row in failure_cases.to_dict("records"):
                                    errors.add_error(row)
                                total_failures += len(failure_cases)
                                exit_code = 1
                                use_logger.error(
                                    "validation_failed",
                                    failures=len(failure_cases),
                                    path=str(failure_path),
                                )

                rows_kept += len(validated_chunk)
                rows_dropped += chunk_rows_total - len(validated_chunk)
                _refresh_column_tracking(validated_chunk)

                if exit_code != 0:
                    aborted = True
                    raise _AbortPipeline(exit_code)

                chunks_emitted = True
                if col_order:
                    yield validated_chunk.reindex(columns=col_order)
                else:
                    yield validated_chunk
        except _AbortPipeline:
            aborted = True
            raise
        except PipelineError:
            raise
        except Exception as exc:  # pragma: no cover - exercised in integration tests
            use_logger.error(
                "chunk_processing_failed",
                error=str(exc),
            )
            aborted = True
            raise _AbortPipeline(1) from exc
        finally:
            if not aborted and not chunks_emitted and not missing_required_columns:
                empty = pd.DataFrame()
                for hook in metadata_hooks:
                    hook_name = _callable_name(hook)
                    try:
                        empty = hook(empty)
                    except Exception as exc:  # pragma: no cover - rare failure path
                        use_logger.error(
                            "metadata_hook_failed",
                            hook=hook_name,
                            error=str(exc),
                            error_type=exc.__class__.__name__,
                            context="empty_frame",
                            strict_mode=strict_mode,
                        )
                        failed_metadata_hooks.add(hook_name)
                        if strict_mode:
                            raise _AbortPipeline(1) from exc
                _refresh_column_tracking(empty)
                yield empty.reindex(columns=col_order) if col_order else empty

    chunk_iterator = _validated_chunks()

    try:
        first_chunk = next(chunk_iterator)
    except StopIteration:
        first_chunk = None
    except _AbortPipeline as abort_exc:
        if missing_required_columns:
            use_logger.warning(
                "validation_skipped",
                missing_columns=sorted(missing_required_columns),
            )
        if total_failures:
            errors.save(failure_path, cfg=cfg)
        else:
            failure_path.unlink(missing_ok=True)
            Path(f"{failure_path}.meta.yaml").unlink(missing_ok=True)
        return abort_exc.code
    except PipelineError:
        if total_failures:
            errors.save(failure_path, cfg=cfg)
        else:
            failure_path.unlink(missing_ok=True)
            Path(f"{failure_path}.meta.yaml").unlink(missing_ok=True)
        return 1

    def _iter_chunks() -> Iterator[pd.DataFrame]:
        if first_chunk is not None:
            yield first_chunk
        yield from chunk_iterator

    csv_path: Path | None = None
    try:
        csv_path = writer(
            _iter_chunks(),
            output_path,
            col_order or None,
            resolved_keys,
        )
    except _AbortPipeline as abort_exc:
        if total_failures:
            errors.save(failure_path, cfg=cfg)
        else:
            failure_path.unlink(missing_ok=True)
            Path(f"{failure_path}.meta.yaml").unlink(missing_ok=True)
        output_path.unlink(missing_ok=True)
        Path(str(output_path) + ".meta.yaml").unlink(missing_ok=True)
        return abort_exc.code
    except Exception as exc:
        if total_failures:
            errors.save(failure_path, cfg=cfg)
        else:
            failure_path.unlink(missing_ok=True)
            Path(f"{failure_path}.meta.yaml").unlink(missing_ok=True)
        output_path.unlink(missing_ok=True)
        Path(str(output_path) + ".meta.yaml").unlink(missing_ok=True)
        use_logger.error(
            "write_fail",
            error=str(exc),
            path=str(output_path),
        )
        return 1

    if total_failures:
        errors.save(failure_path, cfg=cfg)
    else:
        failure_path.unlink(missing_ok=True)
        Path(f"{failure_path}.meta.yaml").unlink(missing_ok=True)

    if optional_cols and (missing_optional := optional_cols - present_columns):
        use_logger.warning(
            "optional_columns_missing",
            columns=sorted(missing_optional),
        )

    if csv_path is None:
        use_logger.error(
            "write_fail", error="writer returned None", path=str(output_path)
        )
        return 1

    resolved_invocation = invocation_tuple

    extra_metadata: dict[str, object] = {}
    if failed_metadata_hooks:
        extra_metadata["metadata_hook_failures"] = sorted(failed_metadata_hooks)

    postprocess_result: PostprocessingPipelineResult | None = None
    if exit_code == 0:
        try:
            postprocess_result = _maybe_run_postprocessing(
                csv_path,
                use_logger,
                postprocess_enabled=effective_postprocess_enabled,
                config_override=effective_postprocess_config,
            )
        except Exception:
            exit_code = 1

    stats_source = stats_extra() if callable(stats_extra) else stats_extra
    stats_payload: dict[str, object] = dict(stats_source or {})

    if postprocess_result is not None:
        postprocess_metadata: dict[str, object] = {
            "output_postprocessed": str(postprocess_result.output_path),
        }
        if postprocess_result.report_path is not None:
            postprocess_metadata["postprocess_report"] = str(
                postprocess_result.report_path
            )
        metrics = postprocess_result.metrics
        if metrics is not None:
            summary = dict(metrics.summary())
            postprocess_metadata["postprocess_metrics"] = summary
            pipeline_version = getattr(metrics, "pipeline_version", None)
            if pipeline_version is not None:
                postprocess_metadata["postprocess_pipeline_version"] = pipeline_version
                stats_payload.setdefault(
                    "postprocess_pipeline_version", pipeline_version
                )
            for key, value in summary.items():
                stats_payload[f"postprocess_{key}"] = value
        extra_metadata.update(postprocess_metadata)

    report = finalise_csv_output(
        csv_path=csv_path,
        rows_total=rows_total,
        rows_kept=rows_kept,
        command=command_str,
        config=config_snapshot,
        inputs=inputs,
        schema=schema_name,
        stats_extra=stats_payload or None,
        invocation=resolved_invocation or None,
        extra_metadata=extra_metadata or None,
    )
    stats = report.stats

    if stats_callback is not None:
        try:
            stats_callback(dict(stats))
        except Exception:  # pragma: no cover - defensive against user callbacks
            use_logger.exception("stats_callback_failed")

    meta_path = report.meta_path

    doc_quality_cfg = getattr(getattr(cfg, "system", None), "doc_quality", None)
    fatal_quality_error = bool(getattr(doc_quality_cfg, "fatal_on_error", False))

    try:
        table_quality(csv_path)
    except Exception as exc:
        tb = traceback.format_exc()
        record_quality_failure(
            meta_path,
            error=str(exc),
            error_type=exc.__class__.__name__,
            traceback=tb,
            fatal=fatal_quality_error,
        )
        log_kwargs = {
            "error": str(exc),
            "error_type": exc.__class__.__name__,
            "path": str(csv_path),
            "traceback": tb,
        }
        if fatal_quality_error:
            log_kwargs["fatal"] = True
            exit_code = 1
        use_logger.warning("quality_report_failed", **log_kwargs)

    if exit_code == 0:
        use_logger.info("write_done", rows=rows_kept, path=str(csv_path))
    return exit_code


def build_parser() -> argparse.ArgumentParser:
    """Create an ``argparse`` parser with shared CSV options.

    Returns
    -------
    argparse.ArgumentParser
        Parser pre-configured with common command-line flags used by
        :mod:`library.utils.cli_tools.csv_utils_main`.
    """
    parser = argparse.ArgumentParser(
        description=(
            "CLI wrapper for :func:`write_csv_deterministic`. This script reads "
            "an input CSV file and re-serialises it deterministically using "
            ":func:`library.common.csv_utils.write_csv_deterministic`."
        )
    )
    add_common_arguments(parser)
    parser.add_argument("--col-order", nargs="*", help="Preferred column order")
    parser.add_argument("--key-cols", nargs="*", help="Columns used for sorting")
    parser.add_argument(
        "--chunk-size",
        type=positive_int,
        default=1000,
        help="Number of rows read per chunk",
    )
    parser.add_argument(
        "--merge-chunk-size",
        type=positive_int,
        default=1000,
        help="Rows loaded per temporary file during merge",
    )
    parser.add_argument(
        "--config",
        dest="config",
        type=path_argument,
        default=DEFAULT_CONFIG_PATH,
        help=f"YAML configuration file (default: {DEFAULT_CONFIG_PATH})",
    )
    parser.add_argument(
        "--print-config",
        action="store_true",
        help="Print effective configuration and exit",
    )
    return parser

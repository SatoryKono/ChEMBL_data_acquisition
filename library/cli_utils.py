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
import uuid
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TypedDict, cast

import numpy as np
import pandas as pd
from pandera.errors import SchemaErrors

from . import cli
from .cli import (
    Logger,
    LoggerConfig,
    add_common_arguments,
    apply_config_overrides,
    path_argument,
    positive_int,
)
from .cli.pipeline_definition import (
    Fetcher,
    MetadataHook,
    PipelineDefinition,
    TableQualityHook,
    Validator,
    Writer,
    normalise_definition,
)
from .cli.run_context import compute_generated_at
from .common.log import logger as default_logger
from .config import Config, ConfigError, ensure_dirs, print_config
from .config.loader import DEFAULT_CONFIG_PATH
from .io import StandardOutputArtifacts, save_standard_outputs
from .maintenance import ensure_legacy_cleanup
from .metadata import Stats, file_sha256, record_quality_failure, write_meta_yaml
from .sidecar import SidecarErrors
from .table_quality import TableQualityProfiler
from .utils.qc_report import build_reports_from_profiler

__all__ = [
    "Fetcher",
    "MetadataHook",
    "PipelineDefinition",
    "PipelineError",
    "PipelineExecutionResult",
    "RunPipelineResult",
    "TableQualityHook",
    "Validator",
    "Writer",
    "build_parser",
    "normalise_definition",
    "resolve_invocation",
    "canonical_run_descriptor",
    "ensure_run_id",
    "run_cli_command",
    "run_pipeline",
]


class StatsExtra(TypedDict, total=False):
    """Additional metrics merged into the metadata output."""

    chunk_fetch_failure_chunks: int
    chunk_fetch_failure_ids: list[str]
    chunk_fetch_failure_ids_total: int
    chunk_fetch_failure_ids_truncated: bool
    custom_metric: int


StatsExtraMapping = Mapping[str, object]


@dataclass(slots=True, frozen=True)
class SchemaSnapshot:
    """Captures introspection details derived from a Pandera schema."""

    required_columns: frozenset[str]
    optional_columns: frozenset[str]
    column_dtypes: dict[str, object | None]

    @classmethod
    def from_schema(cls, schema: object | None) -> SchemaSnapshot:
        if schema is None:
            return cls(
                required_columns=frozenset(),
                optional_columns=frozenset(),
                column_dtypes={},
            )

        schema_columns = getattr(schema, "columns", {}) or {}
        required_cols = {
            name
            for name, column in schema_columns.items()
            if getattr(column, "required", False)
        }
        optional_cols = set(schema_columns) - required_cols
        column_dtypes = {
            name: getattr(column, "dtype", None)
            for name, column in schema_columns.items()
        }
        return cls(
            required_columns=frozenset(required_cols),
            optional_columns=frozenset(optional_cols),
            column_dtypes=column_dtypes,
        )


@dataclass(slots=True)
class PipelineExecutionResult:
    """Outcome of :func:`run_pipeline` including output artefact paths."""

    exit_code: int
    dataset_path: Path | None
    failure_path: Path | None = None
    metadata_path: Path | None = None

    def __int__(self) -> int:  # pragma: no cover - convenience for legacy callers
        return int(self.exit_code)

    def __eq__(self, other: object) -> bool:  # pragma: no cover - legacy helpers
        if isinstance(other, PipelineExecutionResult):
            return self.exit_code == other.exit_code
        if isinstance(other, int):
            return self.exit_code == other
        return NotImplemented


@dataclass(slots=True)
class PipelineMetrics:
    """Aggregates per-chunk counters into totals for metadata."""

    rows_total: int = 0
    rows_kept: int = 0
    rows_dropped: int = 0

    def add_total(self, total: int) -> None:
        self.rows_total += total

    def record(self, *, total: int, kept: int) -> None:
        self.rows_kept += kept
        self.rows_dropped += total - kept

    def build_stats(
        self, *, output_sha256: str, extra: StatsExtraMapping | None
    ) -> Stats:
        stats: Stats = {
            "rows_total": self.rows_total,
            "rows_kept": self.rows_kept,
            "rows_dropped": self.rows_dropped,
            "output_sha256": output_sha256,
        }
        if extra:
            stats.update(extra)
        return stats


@dataclass(slots=True, frozen=True)
class RunPipelineResult:
    """Return value exposing the exit code alongside output artefact paths."""

    # NOTE:
    # Python 3.13 tightened the validation rules around ``__slots__`` when
    # inheriting from built-in types which broke the historical ``int``
    # subclass used by the CLI entry points.  The dataclass-based wrapper keeps
    # the public API intact (``int(result)`` and truthiness checks) without
    # triggering the import error that Windows users observed when running the
    # data acquisition scripts.

    exit_code: int
    dataset_path: Path | None
    artifacts: StandardOutputArtifacts | None = None

    def __int__(self) -> int:
        """Return the stored exit code."""

        return int(self.exit_code)

    def __index__(self) -> int:  # pragma: no cover - mirrors __int__
        return int(self.exit_code)

    def __bool__(self) -> bool:  # pragma: no cover - mirrors ``bool(int)`` semantics
        return bool(self.exit_code)

    def __eq__(self, other: object) -> bool:  # pragma: no cover - trivial delegation
        if isinstance(other, RunPipelineResult):
            return (
                self.exit_code == other.exit_code
                and self.dataset_path == other.dataset_path
                and self.artifacts == other.artifacts
            )
        if isinstance(other, int):
            return self.exit_code == other
        return NotImplemented


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
    path_entries: list[str] = []
    for field in path_fields:
        value = getattr(args, field, None)
        if value in (None, argparse.SUPPRESS):
            continue
        normalised = _normalise_path(value)
        if normalised is not None:
            path_entries.append(f"{field}={normalised}")
    if path_entries:
        parts.extend(sorted(path_entries))
    return "\n".join(parts)


def canonical_run_descriptor(
    args: argparse.Namespace, parser: argparse.ArgumentParser
) -> str:
    """Return a normalised descriptor summarising ``args`` and resolved paths."""

    return _canonical_run_descriptor(args, parser)


def ensure_run_id(
    args: argparse.Namespace, parser: argparse.ArgumentParser, log_cfg: LoggerConfig
) -> str | None:
    """Populate ``args.run_id``/``log_cfg.run_id`` when not explicitly provided."""

    run_id_value = getattr(args, "run_id", None)
    if isinstance(run_id_value, str):
        run_id_value = run_id_value.strip() or None
    elif run_id_value in (argparse.SUPPRESS,):
        run_id_value = None

    if not run_id_value:
        descriptor = canonical_run_descriptor(args, parser)
        if descriptor:
            run_id_value = uuid.uuid5(uuid.NAMESPACE_URL, descriptor).hex
        else:
            candidate = getattr(log_cfg, "run_id", "")
            run_id_value = candidate.strip() or None

    if run_id_value:
        args.run_id = run_id_value
        log_cfg.run_id = run_id_value

    return run_id_value


def run_cli_command(
    *,
    args: argparse.Namespace,
    parser: argparse.ArgumentParser,
    base_parser: argparse.ArgumentParser | None = None,
    log_cfg: LoggerConfig,
    mapping: Mapping[str, str],
    run: Callable[[Config, argparse.Namespace], int],
    logger: Logger | None = None,
) -> int:
    """Execute CLI boilerplate shared by data acquisition commands."""

    debug_flag = bool(getattr(args, "debug", False))
    keep_flag = bool(getattr(args, "keep_intermediate", False))
    legacy_flag = bool(getattr(args, "emit_legacy_artifacts", False))
    diagnostics_enabled = legacy_flag or debug_flag or keep_flag
    if diagnostics_enabled != legacy_flag:
        args.emit_legacy_artifacts = diagnostics_enabled

    desired_level = getattr(args, "log_level", log_cfg.level)
    if getattr(args, "verbose", False):
        desired_level = "DEBUG"
    if not desired_level:
        desired_level = "INFO"
    log_cfg.level = str(desired_level).upper()

    def _configure_logging() -> Logger:
        configured_logger = cli.configure_logger(log_cfg)
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
        if not log_cfg.generated_at:
            seed_parts: list[str] = []
            invocation = getattr(args, "invocation", None)
            if isinstance(invocation, Sequence) and invocation:
                seed_parts.extend(str(part) for part in invocation)
            else:
                program = getattr(parser, "prog", None)
                if program:
                    seed_parts.append(str(program))
                seed_parts.extend(str(part) for part in sys.argv[1:])
            log_cfg.generated_at = compute_generated_at(
                date_token=getattr(args, "date", None),
                run_id=log_cfg.run_id,
                seed_parts=seed_parts,
            )
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
        cli.prepare_io_paths(args)
    except (ValueError, FileNotFoundError) as exc:
        ensure_run_id(args, parser, log_cfg)
        if not log_cfg.generated_at:
            seed_parts = [str(part) for part in getattr(args, "invocation", ())]
            if not seed_parts:
                program = getattr(parser, "prog", None)
                if program:
                    seed_parts.append(str(program))
                seed_parts.extend(str(part) for part in sys.argv[1:])
            log_cfg.generated_at = compute_generated_at(
                date_token=getattr(args, "date", None),
                run_id=log_cfg.run_id,
                seed_parts=seed_parts,
            )
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
    if not log_cfg.generated_at:
        seed_parts: list[str] = []
        invocation = getattr(args, "invocation", None)
        if isinstance(invocation, Sequence) and invocation:
            seed_parts.extend(str(part) for part in invocation)
        else:
            program = getattr(parser, "prog", None)
            if program:
                seed_parts.append(str(program))
            seed_parts.extend(str(part) for part in sys.argv[1:])
        log_cfg.generated_at = compute_generated_at(
            date_token=getattr(args, "date", None),
            run_id=log_cfg.run_id,
            seed_parts=seed_parts,
        )

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
        doc_quality_cfg = getattr(cfg.system, "doc_quality", None)
        if doc_quality_cfg is not None and hasattr(doc_quality_cfg, "enable"):
            doc_quality_cfg.enable = False

    try:
        if getattr(args, "print_config", False):
            print_config(cfg)
            cli.configure_logger(log_cfg)
            use_logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        cleanup_result = ensure_legacy_cleanup(cfg, logger=use_logger)
        if cleanup_result.performed and not cleanup_result.dry_run:
            log_method = (
                use_logger.warning if cleanup_result.removed_count else use_logger.info
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
        if logger is None:
            use_logger = cli.configure_logger(log_cfg)
        else:
            cli.configure_logger(log_cfg)
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

    exit_code = run(cfg, args)
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


def _resolve_standard_output_naming(
    output_path: Path,
    *,
    cfg: Config | None,
    command: str | None,
    invocation: Sequence[str],
) -> tuple[str, str]:
    """Return ``(table_name, date_tag)`` derived from legacy output paths."""

    name = output_path.name
    if name.startswith(".") and name.endswith(".tmp"):
        name = name[1 : -len(".tmp")]
    stem = Path(name).stem if name else output_path.stem
    remainder = stem

    if remainder.startswith("."):
        stripped = remainder.lstrip(".")
        remainder = stripped or remainder

    prefix = "output."
    while remainder.lower().startswith(prefix) and remainder:
        candidate = remainder[len(prefix) :]
        candidate = candidate.lstrip(".")
        if not candidate:
            break
        remainder = candidate

    table_name = remainder
    date_tag: str | None = None
    if "_" in remainder:
        candidate_name, candidate_suffix = remainder.rsplit("_", 1)
        if candidate_suffix.isdigit():
            table_name = candidate_name or table_name
            date_tag = candidate_suffix
        else:
            table_name = remainder
    if not table_name:
        table_name = "dataset"

    if date_tag is None:
        io_cfg = getattr(cfg, "io", None) if cfg is not None else None
        default_prefix = getattr(io_cfg, "default_date_prefix", None)
        seed_parts: list[str] = []
        if command:
            seed_parts.append(command)
        if invocation:
            seed_parts.extend(invocation)
        seed_parts.append(str(output_path))
        generated_at = compute_generated_at(
            date_token=default_prefix,
            run_id=str(output_path),
            seed_parts=seed_parts,
        )
        date_tag = generated_at.split("T", 1)[0].replace("-", "")

    return table_name, date_tag


def run_pipeline(
    *,
    definition: PipelineDefinition | None = None,
    fetcher: Fetcher,
    output_path: Path,
    failure_path: Path,
    cfg: Config | None = None,
    logger: Logger | None = None,
    emit_standard_outputs: bool = True,
    emit_legacy_artifacts: bool = True,
    cleanup_standard_source: bool | None = None,
    **legacy_kwargs: object,
) -> PipelineExecutionResult:
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
    emit_standard_outputs:
        When ``True`` the validated dataset is materialised via
        :func:`library.io.save_standard_outputs` together with QC reports.
    emit_legacy_artifacts:
        When ``True`` the historical deterministic writer and sidecar artefacts are
        produced.
    cleanup_standard_source:
        Optional override controlling whether :func:`library.io.save_standard_outputs`
        removes the original dataset when the canonical path differs. Defaults to
        ``True`` when omitted, matching the historical behaviour.
    legacy_kwargs:
        Deprecated keyword arguments matching the historical signature. When
        provided they are converted into a :class:`PipelineDefinition`.

    Returns
    -------
    PipelineExecutionResult
        Object containing the exit code together with the resolved dataset
        path and optional standard output artifacts. Zero indicates success while
        non-zero codes signal validation or writer failures.
        ``PipelineError`` exceptions raised by ``fetcher`` are converted into a non-zero exit code.
    """

    use_logger = logger or default_logger

    standard_outputs_enabled = bool(emit_standard_outputs)
    legacy_outputs_enabled = bool(emit_legacy_artifacts)

    definition = normalise_definition(definition, legacy_kwargs)

    schema = definition.schema
    schema_name = definition.schema_name
    metadata_hooks = list(definition.metadata_hooks)
    validators = list(definition.validators)
    writer = definition.writer
    config_snapshot = dict(definition.config_snapshot)
    inputs = dict(definition.inputs)
    key_columns = list(definition.key_columns)
    table_quality = definition.table_quality or (lambda _: None)
    stats_extra_config = definition.stats_extra
    stats_callback = definition.stats_callback
    strict_mode = bool(definition.strict_mode)
    invocation_tuple = tuple(str(part) for part in (definition.invocation or ()))
    command = definition.command
    dictionary_resources = tuple(definition.dictionary_resources or ()) or None

    profiler: TableQualityProfiler | None = (
        TableQualityProfiler() if standard_outputs_enabled else None
    )
    collected_chunks: list[pd.DataFrame] | None = (
        [] if standard_outputs_enabled else None
    )
    standard_artifacts: StandardOutputArtifacts | None = None
    dataset_path: Path | None = None

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

    schema_snapshot = SchemaSnapshot.from_schema(schema)
    required_cols = set(schema_snapshot.required_columns)
    optional_cols = set(schema_snapshot.optional_columns)
    schema_column_dtypes = dict(schema_snapshot.column_dtypes)

    errors = SidecarErrors()
    meta_path: Path | None = None
    metrics = PipelineMetrics()
    total_failures = 0
    exit_code = 0
    present_columns: set[str] = set()
    missing_required_columns: set[str] = set()
    all_columns: set[str] = set()
    validation_enabled = True
    failed_metadata_hooks: set[str] = set()

    def _cleanup_failures(has_failures: bool) -> None:
        if emit_legacy_artifacts and has_failures:
            errors.save(failure_path, cfg=cfg)
        else:
            failure_path.unlink(missing_ok=True)
            Path(f"{failure_path}.meta.yaml").unlink(missing_ok=True)

    def _build_result(
        code: int, *, dataset: Path | None = None
    ) -> PipelineExecutionResult:
        failure_artifact = failure_path if emit_legacy_artifacts else None
        metadata_artifact = meta_path if emit_legacy_artifacts else None
        return PipelineExecutionResult(
            exit_code=code,
            dataset_path=dataset,
            failure_path=failure_artifact,
            metadata_path=metadata_artifact,
        )

    try:
        iterable = _as_iterable(fetcher())
    except PipelineError:
        _cleanup_failures(False)
        return _build_result(1)
    except Exception as exc:  # pragma: no cover - exercised in integration tests
        use_logger.error("fetch_failed", error=str(exc), exc_info=exc)
        _cleanup_failures(False)
        return _build_result(1)

    class _AbortPipeline(RuntimeError):
        """Internal exception raised to abort processing early."""

        def __init__(self, code: int = 1) -> None:
            super().__init__("pipeline aborted")
            self.code = code

    schema_columns: list[str] | None
    optional_column_order: list[str]
    if schema is not None:
        schema_columns = list(schema_column_dtypes)
        optional_column_order = [
            column for column in schema_columns if column in optional_cols
        ]
    else:
        schema_columns = None
        optional_column_order = []

    col_order: list[str] = list(schema_columns or [])
    resolved_keys: list[str] = []

    def _refresh_column_tracking(frame: pd.DataFrame) -> None:
        chunk_columns = set(frame.columns)
        present_columns.update(chunk_columns)
        all_columns.update(chunk_columns)

        if optional_column_order:
            all_columns.update(optional_column_order)

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
        nonlocal metrics, total_failures, exit_code
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
                metrics.add_total(chunk_rows_total)

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
                            exc_info=exc,
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
                                    exc_info=exc,
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

                metrics.record(total=chunk_rows_total, kept=len(validated_chunk))
                _refresh_column_tracking(validated_chunk)

                if exit_code != 0:
                    aborted = True
                    raise _AbortPipeline(exit_code)

                chunks_emitted = True
                source_columns = list(validated_chunk.columns)
                if col_order:
                    emitted_chunk = validated_chunk.reindex(columns=col_order)
                else:
                    emitted_chunk = validated_chunk
                emitted_chunk.attrs["_source_columns"] = source_columns
                yield emitted_chunk
        except _AbortPipeline:
            aborted = True
            raise
        except PipelineError:
            raise
        except Exception as exc:  # pragma: no cover - exercised in integration tests
            use_logger.error(
                "chunk_processing_failed",
                error=str(exc),
                exc_info=exc,
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
                            exc_info=exc,
                        )
                        failed_metadata_hooks.add(hook_name)
                        if strict_mode:
                            raise _AbortPipeline(1) from exc
                _refresh_column_tracking(empty)
                emitted_empty = empty.reindex(columns=col_order) if col_order else empty
                emitted_empty.attrs["_source_columns"] = list(empty.columns)
                yield emitted_empty

    optional_cols_added: set[str] = set()

    def _resolve_dtype(dtype: object | None) -> object | None:
        if dtype is None:
            return None
        dtype_type = getattr(dtype, "type", None)
        if dtype_type is not None:
            try:
                return np.dtype(dtype_type)
            except (TypeError, ValueError):
                pass
        python_type = getattr(dtype, "python_type", None)
        if isinstance(python_type, type):
            return python_type
        if isinstance(dtype, type) or isinstance(dtype, str):
            return dtype
        return dtype

    def _prepare_chunk(frame: pd.DataFrame) -> pd.DataFrame:
        nonlocal optional_cols_added

        prepared = frame
        original_columns = list(frame.attrs.get("_source_columns", frame.columns))
        original_set = set(original_columns)
        if optional_column_order:
            missing = [
                column for column in optional_column_order if column not in original_set
            ]
            if missing:
                optional_cols_added.update(missing)
                existing = list(prepared.columns)
                missing_new = [
                    column for column in missing if column not in prepared.columns
                ]
                if missing_new:
                    prepared = prepared.reindex(
                        columns=existing + missing_new, fill_value=""
                    )
                    existing = list(prepared.columns)
                fill_columns = [
                    column for column in missing if column in prepared.columns
                ]
                if fill_columns:
                    fills = {
                        column: prepared[column].fillna("") for column in fill_columns
                    }
                    prepared = prepared.assign(**fills)

        if schema_column_dtypes:
            conversions: dict[str, object] = {}
            for column, dtype in schema_column_dtypes.items():
                if column not in prepared.columns:
                    continue
                resolved_dtype = _resolve_dtype(dtype)
                if resolved_dtype is None:
                    continue
                conversions[column] = resolved_dtype
            if conversions:
                prepared = prepared.astype(conversions)

        prepared.attrs["_source_columns"] = original_columns

        return prepared

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
        _cleanup_failures(total_failures > 0)
        return _build_result(abort_exc.code)
    except PipelineError:
        _cleanup_failures(total_failures > 0)
        return _build_result(1)

    prepared_first: pd.DataFrame | None = None
    if first_chunk is not None:
        prepared_first = _prepare_chunk(first_chunk)
        if profiler is not None:
            profiler.consume(prepared_first)
        if collected_chunks is not None:
            collected_chunks.append(prepared_first)

    def _iter_chunks() -> Iterator[pd.DataFrame]:
        if prepared_first is not None:
            yield prepared_first
        for chunk in chunk_iterator:
            prepared = _prepare_chunk(chunk)
            if profiler is not None:
                profiler.consume(prepared)
            if collected_chunks is not None:
                collected_chunks.append(prepared)
            yield prepared

    csv_path: Path | None = None
    chunk_stream = _iter_chunks()
    should_invoke_writer = legacy_outputs_enabled or not standard_outputs_enabled

    try:
        if should_invoke_writer:
            csv_path = writer(
                chunk_stream,
                output_path,
                col_order or None,
                resolved_keys,
            )
        else:
            for _ in chunk_stream:
                pass
    except _AbortPipeline as abort_exc:
        output_path.unlink(missing_ok=True)
        Path(str(output_path) + ".meta.yaml").unlink(missing_ok=True)
        _cleanup_failures(total_failures > 0)
        return _build_result(abort_exc.code)
    except Exception as exc:
        output_path.unlink(missing_ok=True)
        Path(str(output_path) + ".meta.yaml").unlink(missing_ok=True)
        use_logger.error(
            "write_fail",
            error=str(exc),
            path=str(output_path),
            exc_info=exc,
        )
        _cleanup_failures(total_failures > 0)
        return _build_result(1)

    if optional_column_order:
        present_columns.update(optional_cols_added)

    _cleanup_failures(total_failures > 0)

    if optional_cols and (missing_optional := optional_cols - present_columns):
        use_logger.warning(
            "optional_columns_missing",
            columns=sorted(missing_optional),
        )

    # Remove legacy artifacts' meta file if not emitting legacy artifacts
    if not emit_legacy_artifacts:
        Path(str(output_path) + ".meta.yaml").unlink(missing_ok=True)

    # Standard outputs logic
    if standard_outputs_enabled and exit_code == 0:
        frames_to_merge = list(collected_chunks or [])
        if frames_to_merge:
            dataset_frame = pd.concat(frames_to_merge, ignore_index=True)
        else:
            empty_columns = list(col_order) if col_order else sorted(all_columns)
            dataset_frame = pd.DataFrame(columns=empty_columns)
        if col_order:
            dataset_frame = dataset_frame.reindex(columns=col_order)
        elif all_columns and list(dataset_frame.columns) != sorted(all_columns):
            dataset_frame = dataset_frame.reindex(columns=sorted(all_columns))

        active_profiler = profiler or TableQualityProfiler()
        quality_report, correlation_report = build_reports_from_profiler(
            active_profiler
        )

        table_name, date_tag = _resolve_standard_output_naming(
            output_path,
            cfg=cfg,
            command=command_str,
            invocation=invocation_tuple,
        )

        io_cfg = getattr(cfg, "io", None) if cfg is not None else None
        output_dir_value = (
            getattr(io_cfg, "output_dir", None) if io_cfg is not None else None
        )
        if output_dir_value is None:
            raise ValueError(
                "cfg.io.output_dir must be configured for standard outputs"
            )

        cleanup_canonical_source = (
            bool(cleanup_standard_source)
            if cleanup_standard_source is not None
            else True
        )

        standard_artifacts = save_standard_outputs(
            dataset_frame,
            correlation_report,
            quality_report,
            table_name=table_name,
            date_tag=date_tag,
            output_path=output_path,
            output_dir=Path(output_dir_value),
            key_columns=key_columns,
            cleanup_source=cleanup_canonical_source,
        )
        use_logger.info(
            "standard_outputs_written",
            dataset=str(standard_artifacts.dataset),
            quality_report=str(standard_artifacts.quality_report),
            correlation_report=str(standard_artifacts.correlation_report),
        )
        dataset_path = standard_artifacts.dataset
        csv_path = standard_artifacts.dataset

    if dataset_path is None and csv_path is not None:
        dataset_path = csv_path

    quality_input_path: Path | None
    if dataset_path is not None:
        quality_input_path = Path(dataset_path)
    elif csv_path is not None:
        quality_input_path = Path(csv_path)
    else:
        quality_input_path = None

    if legacy_outputs_enabled and csv_path is None:
        use_logger.error(
            "write_fail", error="writer returned None", path=str(output_path)
        )
        return _build_result(1)

    if dataset_path is None:
        use_logger.error(
            "write_fail", error="no dataset produced", path=str(output_path)
        )
        return _build_result(1)

    if stats_extra_config is None:
        extra_stats: StatsExtraMapping | None = None
    elif callable(stats_extra_config):
        extra_candidate = cast(
            StatsExtra | StatsExtraMapping | None, stats_extra_config()
        )
        extra_stats = dict(extra_candidate) if extra_candidate is not None else None
    else:
        extra_stats = dict(cast(StatsExtra | StatsExtraMapping, stats_extra_config))

    stats = metrics.build_stats(
        output_sha256=file_sha256(dataset_path),
        extra=extra_stats,
    )
    stats["dataset_path"] = str(dataset_path)

    if stats_callback is not None:
        try:
            stats_callback(dict(stats))
        except Exception:  # pragma: no cover - defensive against user callbacks
            use_logger.exception("stats_callback_failed")

    resolved_invocation = invocation_tuple

    meta_path: Path | None = None
    extra_metadata: dict[str, object] | None = None
    if legacy_outputs_enabled and failed_metadata_hooks:
        extra_metadata = {"metadata_hook_failures": sorted(failed_metadata_hooks)}

    if emit_legacy_artifacts and quality_input_path is not None:
        meta_path = write_meta_yaml(
            csv_path=quality_input_path,
            command=command_str,
            config=config_snapshot,
            inputs=inputs,
            stats=stats,
            schema=schema_name,
            invocation=resolved_invocation or None,
            extra_metadata=extra_metadata or None,
            dictionary_resources=dictionary_resources,
        )
        doc_quality_cfg = getattr(getattr(cfg, "system", None), "doc_quality", None)
        fatal_quality_error = bool(getattr(doc_quality_cfg, "fatal_on_error", False))

        quality_source = dataset_path or csv_path
        try:
            if quality_source is not None:
                table_quality(quality_source)
        except Exception as exc:
            tb = traceback.format_exc()
            record_quality_failure(
                meta_path,
                error=str(exc),
                error_type=exc.__class__.__name__,
                traceback=tb,
                fatal=fatal_quality_error,
            )
            fallback_path = quality_source if quality_source is not None else csv_path
            if fallback_path is None:
                fallback_path = output_path
            log_kwargs = {
                "error": str(exc),
                "error_type": exc.__class__.__name__,
                "path": str(fallback_path),
                "traceback": tb,
            }
            if fatal_quality_error:
                log_kwargs["fatal"] = True
                exit_code = 1
            use_logger.warning("quality_report_failed", **log_kwargs)

    if exit_code == 0:
        write_path = dataset_path or csv_path or output_path
        use_logger.info("write_done", rows=metrics.rows_kept, path=str(write_path))
    return RunPipelineResult(exit_code, dataset_path, standard_artifacts)


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

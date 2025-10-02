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
import sys
import traceback
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from pathlib import Path
from typing import Protocol, TypeVar, overload

import shlex

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
from .config import Config, ensure_dirs, print_config
from .log import logger as default_logger
from .metadata import Stats, file_sha256, write_meta_yaml
from .sidecar import SidecarErrors
from .utils.config import DEFAULT_CONFIG_PATH

SchemaT = TypeVar("SchemaT")

class ValidationResult(Protocol):
    """Protocol describing the return type of validator callables."""

    data: pd.DataFrame
    failure_cases: pd.DataFrame


class Validator(Protocol):
    """Callable interface for data frame validators."""

    def __call__(
        self, df: pd.DataFrame
    ) -> ValidationResult:  # pragma: no cover - Protocol
        ...


MetadataHook = Callable[[pd.DataFrame], pd.DataFrame]
Writer = Callable[
    [Iterable[pd.DataFrame], Path, Sequence[str] | None, Sequence[str]], Path
]
TableQualityHook = Callable[[Path], None]
Fetcher = Callable[[], Iterable[pd.DataFrame] | pd.DataFrame]


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

    log_cfg.level = getattr(args, "log_level", log_cfg.level)
    configured_logger = cli.configure_logger(log_cfg)
    use_logger = logger or configured_logger
    use_logger.info("pipeline_start", run_id=log_cfg.run_id)

    try:
        config_arg = getattr(args, "config", None)
        if isinstance(config_arg, (str, Path)):
            config_path = config_arg
        else:
            default_config = parser.get_default("config")
            if not isinstance(default_config, (str, Path)):
                msg = "configuration path must be provided"
                raise ValueError(msg)
            config_path = default_config

        cli.prepare_io_paths(args)

        cfg: Config = apply_config_overrides(
            args,
            parser,
            config_path,
            mapping=dict(mapping),
            base_parser=base_parser,
        )
        if getattr(args, "print_config", False):
            print_config(cfg)
            cli.configure_logger(log_cfg)
            use_logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        use_logger = cli.configure_logger(log_cfg)
    except (ValueError, TypeError) as exc:
        use_logger.error(
            "config_error",
            error=str(exc),
            config=str(getattr(args, "config", "")),
        )
        use_logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        use_logger.error(
            "directory_setup_failed",
            error=str(exc),
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


def run_pipeline(
    *,
    fetcher: Fetcher,
    schema: SchemaT | None,
    schema_name: str,
    validators: Sequence[Validator] | None,
    metadata_hooks: Sequence[MetadataHook] | None,
    writer: Writer,
    output_path: Path,
    failure_path: Path,
    command: str | None = None,
    invocation: Sequence[str] | None = None,
    config_snapshot: Mapping[str, object],
    inputs: Mapping[str, object],
    key_columns: Sequence[str],
    table_quality: TableQualityHook,
    cfg: Config | None = None,
    strict_mode: bool = False,
    logger: Logger | None = None,
) -> int:
    """Execute a data pipeline and write deterministic CSV output.

    Parameters
    ----------
    fetcher:
        Callable returning an iterable of raw :class:`pandas.DataFrame`
        objects.  Each frame represents a chunk of data retrieved from an
        upstream service.
    schema:
        Pandera schema describing the expected output columns.  The helper
        inspects ``schema.columns`` to determine required and optional fields
        as well as to construct the preferred column order.
    schema_name:
        Human readable name of ``schema`` persisted in the metadata file.
    validators:
        Sequence of callables returning Pandera ``ValidationResult`` objects.
        Every validated chunk is passed through each validator in order.
    metadata_hooks:
        Callables applied sequentially to every chunk before validation.
    writer:
        Function responsible for serialising the validated chunks to ``CSV``.
        It receives the chunk iterator, destination path, final column order
        and the subset of ``key_columns`` present in the output.
    output_path:
        Destination ``CSV`` file.
    failure_path:
        Path for persisting validation failure cases.
    command:
        Command used to launch the pipeline.  Stored in metadata output.
    invocation:
        Optional command invocation captured as a sequence of arguments. When
        provided it is persisted to metadata alongside the joined ``command``
        string.
    config_snapshot:
        Mapping of configuration values persisted to metadata.
    inputs:
        Mapping of input file descriptions persisted to metadata.
    key_columns:
        Preferred key columns used when sorting deterministic output.  Only
        columns present in the final dataset are forwarded to ``writer``.
    table_quality:
        Callable invoked after writing the CSV to compute quality metrics.
    cfg:
        Optional application configuration forwarded to sidecar metadata.
    strict_mode:
        When ``True`` metadata hook failures abort the pipeline mimicking the
        legacy behaviour used in CI environments.
    logger:
        Optional logger.  Defaults to :data:`library.log.logger` when omitted.

    Returns
    -------
    int
        Zero on success, one when validation fails or ``writer`` raises an
        exception.  ``PipelineError`` exceptions raised by ``fetcher`` are
        converted into a non-zero return code.
    """

    use_logger = logger or default_logger

    # NOTE:
    # ``invocation`` is an optional parameter which, in practice, might be
    # omitted by older call-sites.  When that happens Python still provides the
    # default ``None`` value, however certain execution environments (for
    # example when the function is referenced through ``functools.partial`` or
    # dynamically inspected) have been observed to trigger ``NameError`` while
    # the default is being resolved.  To keep the metadata handling resilient we
    # retrieve the argument from ``locals()`` instead of referencing the name
    # directly which guarantees the lookup succeeds even when the optimiser
    # elides the symbol.
    invocation_value = locals().get("invocation")

    invocation_tuple: tuple[str, ...] = ()
    if invocation_value is not None:
        invocation_tuple = tuple(str(part) for part in invocation_value)

    if command is not None:
        command_str = command
    elif invocation_tuple:
        command_str = " ".join(shlex.quote(part) for part in invocation_tuple)
    else:
        raise ValueError("run_pipeline requires either 'command' or 'invocation'")

    metadata_hooks = list(metadata_hooks or [])
    validators = list(validators or [])

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
        use_logger.error("write_fail", error="writer returned None", path=str(output_path))
        return 1

    stats: Stats = {
        "rows_total": rows_total,
        "rows_kept": rows_kept,
        "rows_dropped": rows_dropped,
        "output_sha256": file_sha256(csv_path),
    }
 
    resolved_invocation = invocation_tuple

    extra_metadata: dict[str, object] = {}
    if failed_metadata_hooks:
        extra_metadata["metadata_hook_failures"] = sorted(failed_metadata_hooks)

    meta_path = write_meta_yaml(
        csv_path=csv_path,
        command=command_str,
        config_subset=config_snapshot,
        inputs=inputs,
        stats=stats,
        schema=schema_name,
        invocation=resolved_invocation or None,
        extra_metadata=extra_metadata or None,
    )

    try:
        table_quality(csv_path)
    except Exception as exc:
        use_logger.error(
            "quality_report_failed",
            error=str(exc),
            error_type=exc.__class__.__name__,
            path=str(csv_path),
            traceback=traceback.format_exc(),

        )
        if schema is not None:
            Path(csv_path).unlink(missing_ok=True)
        meta_path.unlink(missing_ok=True)
        return 1

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
            ":func:`library.csv_utils.write_csv_deterministic`."
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

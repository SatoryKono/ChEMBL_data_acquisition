"""Command line interface for retrieving ChEMBL assay data.

The script provides a reusable :func:`main` entry point together with helper
functions that unit tests import directly. Functions return explicit exit codes
instead of terminating the interpreter to make orchestration easier.
"""

from __future__ import annotations

import argparse
import sys
from collections import deque
from collections.abc import Iterable, Iterator, Mapping, Sequence
from functools import partial
from itertools import islice
from pathlib import Path
from time import sleep
from typing import Any

import pandas as pd
import requests

from library import cli, io
from library.cli import (
    ConfigMetadata,
    Logger,
    LoggerConfig,
    configure_logger,
)
from library.cli import build_parser as base_parser
from library.cli.base import PipelineCLIBase
from library.cli.metadata import prepare_option
from library.cli.pipeline_definition import PipelineDefinition
from library.cli_utils import run_cli_command, run_pipeline
from library.common.csv_utils import write_csv_chunks_deterministic
from library.common.fetch_retry import ChunkFailureTracker, compute_backoff_delay
from library.common.log import logger
from library.config import Config, _serialize_paths
from library.integration import chembl_library as cl
from library.orchestration import ETLContext
from library.pipelines.assay import AssayPipelineOptions
from library.pipelines.assay.chembl_assay import ASSAY_COLUMNS, MAX_ASSAY_CHUNK_SIZE
from library.pipelines.common import (
    ChunkedFetchConfig,
    CsvWriterConfig,
    PipelineRunResult,
    add_pipeline_metadata,
    prepare_chunked_pipeline,
)
from library.pipelines.common.metadata import get_pipeline_version
from library.qa.reporting import build_table_quality_hook
from library.utils.data_correlation import generate_correlation_report
from library.utils.qc_report import generate_qc_report
from library.schemas import AssaysSchema, normalize_assays
from library.validation import validate_assays as validate_assays_schema

__all__ = [
    "configure_logger",
    "main",
    "run_assay_service",
    "run",
    "run_chembl",
    "run_cli_command",
]


DEFAULT_INPUT_NAME = "assay.csv"
DEFAULT_OUTPUT_STEM = "assays"

# Backwards compatibility: legacy configs referenced the private
# ``_ASSAY_MAX_IDS_PER_REQUEST`` constant before it was renamed to
# :data:`MAX_ASSAY_CHUNK_SIZE`.  Re-expose the alias so that pipelines relying on
# the old setting fail gracefully instead of raising ``NameError`` during chunk
# processing.
ASSAY_MAX_IDS_PER_REQUEST = MAX_ASSAY_CHUNK_SIZE
_ASSAY_MAX_IDS_PER_REQUEST = MAX_ASSAY_CHUNK_SIZE

_ASSAY_OUTPUT_DROP_COLUMNS = [
    "ASSAY_ID",
    "Target TYPE",
    "acts_per_assay_step5",
    "cited_assay_corr",
    "error_assay_corr",
    "higly_correlated_cit",
    "month",
    "shuffled_cit",
    "shuffled_target_assay",
    "substrate_name",
    "target_name",
    "version",
    "assay_category",
    "src_assay_id",
]


def _drop_assay_output_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Remove disallowed columns from the final assay output."""

    allowed_cols = [
        column for column in frame.columns if column not in _ASSAY_OUTPUT_DROP_COLUMNS
    ]
    dropped_present = [
        column for column in _ASSAY_OUTPUT_DROP_COLUMNS if column in frame.columns
    ]

    # ``errors='ignore'`` guarantees the pipeline remains stable if a column is absent.
    trimmed = frame.drop(columns=_ASSAY_OUTPUT_DROP_COLUMNS, errors="ignore")
    trimmed = trimmed.loc[:, allowed_cols]

    dropped_display = ", ".join(dropped_present) if dropped_present else "<none>"
    logger.log(
        "INFO",
        "get_assay_data",
        msg=f"Dropped columns from output.assay_*: {dropped_display}",
        dropped_columns=dropped_present,
    )

    return trimmed


ASSAY_OUTPUT_DROP_COLUMNS: list[str] = [
    "ASSAY_ID",
    "Target TYPE",
    "acts_per_assay_step5",
    "cited_assay_corr",
    "error_assay_corr",
    "higly_correlated_cit",
    "month",
    "shuffled_cit",
    "shuffled_target_assay",
    "substrate_name",
    "target_name",
    "version",
    "assay_category",
    "src_assay_id",
]


def remove_assay_output_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Return ``df`` without columns disallowed in ``output.assay_*`` exports."""

    allowed_cols = [
        column for column in df.columns if column not in ASSAY_OUTPUT_DROP_COLUMNS
    ]
    cleaned = df.drop(columns=ASSAY_OUTPUT_DROP_COLUMNS, errors="ignore")
    if allowed_cols:
        cleaned = cleaned.loc[:, allowed_cols]
    return cleaned


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute assay retrieval from the ChEMBL API.

    The output CSV arranges columns so that fields defined in
    :class:`~library.schemas.assays.AssaysSchema` appear first. Any additional columns
    are appended alphabetically.

    Parameters
    ----------
    cfg : Config
        Application configuration providing API credentials, retry strategy and
        CSV export settings.
    args : argparse.Namespace
        Parsed command-line arguments accepted by the ``assay`` CLI command.

    Returns
    -------
    int
        ``0`` on success, non-zero when validation or I/O failures occur.
        Upstream API errors are logged and converted into a failure code by
        :func:`library.cli_utils.run_pipeline`.
    """
    limit = cfg.assay.limit
    if limit is not None and limit < 0:
        logger.error("invalid_limit", section="assay.limit", limit=limit)
        return 1

    try:
        ids_iter = io.read_ids(args.input_csv, column=cfg.assay.column, cfg=cfg.io)
    except (FileNotFoundError, ValueError) as exc:
        logger.error(
            "read_fail",
            error=str(exc),
            exc_info=exc,
            path=str(args.input_csv),
        )
        return 1

    offset = getattr(args, "offset", 0)
    final_out_attr = getattr(args, "final_out", None)
    if final_out_attr in (None, argparse.SUPPRESS):
        legacy_output = getattr(args, "output_csv", None)
        if legacy_output not in (None, argparse.SUPPRESS):
            output_path = Path(legacy_output)
            if not isinstance(legacy_output, Path):
                args.final_out = output_path
            args.output_csv = output_path
        else:
            output_path = Path(
                io.default_output_path(
                    args.input_csv,
                    cfg.io,
                    date=getattr(args, "date", None),
                )
            )
            args.final_out = output_path
            args.output_csv = output_path
    else:
        output_path = Path(final_out_attr)
        if not isinstance(final_out_attr, Path):
            args.final_out = output_path
        args.output_csv = output_path
    metadata_obj = getattr(args, "_config_metadata", None)
    if not isinstance(metadata_obj, ConfigMetadata):
        metadata_obj = None
    output_source = "cli" if getattr(args, "final_out", None) else "derived"
    logger.info(
        "assay_pipeline_start",
        input=prepare_option(
            metadata_obj, value=str(args.input_csv), default_source="cli"
        ),
        output=prepare_option(
            metadata_obj,
            value=str(output_path),
            default_source=output_source,
        ),
        limit=prepare_option(
            metadata_obj,
            argument="limit",
            path="sources.chembl.pipelines.assay.limit",
            value=limit,
        ),
        offset=prepare_option(
            metadata_obj,
            argument="offset",
            path="sources.chembl.pipelines.assay.offset",
            value=offset,
        ),
        batch_size=prepare_option(
            metadata_obj,
            argument="batch_size",
            path="sources.chembl.pipelines.assay.batch_size",
            value=cfg.assay.batch_size,
        ),
        timeout=prepare_option(
            metadata_obj,
            argument="timeout",
            path="sources.chembl.pipelines.assay.timeout",
            value=cfg.assay.timeout,
        ),
    )
    if offset:
        ids_iter = islice(ids_iter, offset, None)
        logger.info("process_offset", offset=offset)

    processed_ids = 0

    def _iter_ids() -> Iterator[str]:
        nonlocal processed_ids
        for identifier in ids_iter:
            processed_ids += 1
            yield identifier

    if limit is not None:
        ids_source: Iterable[str] = islice(_iter_ids(), limit)
    else:
        ids_source = _iter_ids()

    failure_path = output_path.with_name(f"{output_path.stem}_failure_cases.csv")
    fetch_failure_path = output_path.with_name(f"{output_path.stem}_fetch_failures.csv")

    dropped_columns_seen: set[str] = set()

    def _drop_output_columns(frame: pd.DataFrame) -> pd.DataFrame:
        removed = [
            column for column in ASSAY_OUTPUT_DROP_COLUMNS if column in frame.columns
        ]
        if removed:
            dropped_columns_seen.update(removed)
        return remove_assay_output_columns(frame)

    postprocess_enabled = bool(getattr(args, "postprocess", False))
    emit_legacy = bool(getattr(args, "emit_legacy_artifacts", False))
    dictionary_resources: tuple[str, ...] | None = None

    metadata_hooks = [
        normalize_assays,
        add_pipeline_metadata,
        _drop_output_columns,
        _drop_assay_output_columns,
    ]

    validators: list[Any] = []

    postprocess_runtime = None
    postprocess_runner = None
    postprocess_csv_cfg = None
    postprocess_generate_metrics = None
    postprocess_pipeline_config = None
    SchemaValidationError = None
    StepError = None

    if postprocess_enabled:
        from library.postprocessing import enrich_assay_metadata
        from library.pipelines.assay import postprocessing as ap  # noqa: WPS433
        from library.postprocessing.assays import (
            ASSAY_SCHEMA as POSTPROCESS_ASSAY_SCHEMA,
            run_assay_pipeline as run_assay_postprocess,
            validate_assays as validate_postprocess_assays,
        )
        from library.postprocessing.common.types import (  # noqa: WPS433
            SchemaValidationError as PostprocessSchemaError,
            StepError as PostprocessStepError,
        )
        from library.postprocess import (  # noqa: WPS433
            PostprocessingPipelineConfig,
            generate_metrics_report,
            get_csv_runtime_config,
            get_pipeline_config,
            run_postprocessing_pipeline,
        )

        dictionary_resources = ("dictionary_root",)

        def _enrich_with_dictionaries(frame: pd.DataFrame) -> pd.DataFrame:
            return enrich_assay_metadata(
                frame, dictionary_dir=cfg.resources.dictionary_dir
            )

        metadata_hooks = [
            _enrich_with_dictionaries,
            ap.postprocess_assays,
            *metadata_hooks,
        ]

        validators.append(partial(validate_assays_schema, return_result=True))

        postprocess_pipeline_config = get_pipeline_config(
            "assays", getattr(args, "config", None)
        )
        postprocess_csv_cfg = get_csv_runtime_config(postprocess_pipeline_config)
        postprocess_runtime = PostprocessingPipelineConfig(
            pipeline_config=postprocess_pipeline_config,
            csv_runtime_config=postprocess_csv_cfg,
            runner=run_assay_postprocess,
            validator=validate_postprocess_assays,
            schema=POSTPROCESS_ASSAY_SCHEMA,
            logger=logger,
        )
        postprocess_runner = run_assay_postprocess
        postprocess_generate_metrics = generate_metrics_report
        SchemaValidationError = PostprocessSchemaError
        StepError = PostprocessStepError

    doc_quality_cfg = cfg.system.doc_quality
    if emit_legacy:
        table_quality = build_table_quality_hook(
            doc_quality_cfg,
            table_name=output_path.with_suffix(""),
            destination=output_path.parent,
        )
    else:
        table_quality = lambda _: None  # type: ignore[assignment]

    def _persist_standard_outputs(dataset_csv: Path) -> io.StandardOutputArtifacts:
        dataset_frame = pd.read_csv(
            dataset_csv,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
        )
        table_label = dataset_csv.with_suffix("")
        try:
            correlation_report = generate_correlation_report(
                dataset_frame,
                table_name=str(table_label),
            )
        except Exception as exc:  # pragma: no cover - defensive logging
            logger.warning(
                "assay_correlation_generation_failed",
                error=str(exc),
                path=str(dataset_csv),
            )
            correlation_report = pd.DataFrame()
        try:
            quality_report = generate_qc_report(
                dataset_frame,
                table_name=str(table_label),
            )
        except Exception as exc:  # pragma: no cover - defensive logging
            logger.warning(
                "assay_quality_generation_failed",
                error=str(exc),
                path=str(dataset_csv),
            )
            quality_report = pd.DataFrame()
        table_name_value, date_tag = io.derive_output_labels(
            dataset_csv,
            default_table=DEFAULT_OUTPUT_STEM,
        )
        artifacts = io.save_standard_outputs(
            dataset_frame,
            correlation_report,
            quality_report,
            table_name=table_name_value,
            date_tag=date_tag,
            output_dir=Path(cfg.io.output_dir),
        )
        logger.info(
            "assay_standard_outputs",
            dataset=str(artifacts.dataset),
            quality_report=str(artifacts.quality_report),
            correlation_report=str(artifacts.correlation_report),
        )
        return artifacts

    with ETLContext(cfg) as context:
        client = context.chembl_client

        retry_cfg = cfg.retry
        jitter = retry_cfg.build_jitter()
        chunk_failures = ChunkFailureTracker()

        def fetch_chunk(chunk_ids: Sequence[str]) -> pd.DataFrame:
            attempts = max(1, retry_cfg.max_attempts)
            pending: deque[list[str]] = deque([list(chunk_ids)])
            frames: list[pd.DataFrame] = []

            while pending:
                current = pending.popleft()
                if not current:
                    continue

                for attempt in range(1, attempts + 1):
                    try:
                        frame = cl.get_assays(
                            current,
                            cfg=cfg.api,
                            client=client,
                            chunk_size=min(cfg.assay.batch_size, len(current)),
                            timeout=cfg.assay.timeout,
                        )
                    except (requests.RequestException, ValueError) as exc:
                        error_message = str(exc)
                        context = {
                            "chunk_ids": list(current),
                            "chunk_size": len(current),
                            "attempt": attempt,
                            "max_attempts": attempts,
                            "batch_size": cfg.assay.batch_size,
                            "timeout": cfg.assay.timeout,
                        }
                        log_context = {
                            k: v for k, v in context.items() if k != "chunk_ids"
                        }

                        if attempt >= attempts:
                            if len(current) > 1:
                                split_index = max(1, len(current) // 2)
                                left = current[:split_index]
                                right = current[split_index:]
                                logger.warning(
                                    "assay_fetch_split",
                                    extra={
                                        "msg": error_message,
                                        "chunk_ids": context["chunk_ids"],
                                    },
                                    **log_context,
                                )
                                if right:
                                    pending.appendleft(right)
                                if left:
                                    pending.appendleft(left)
                                break

                            logger.error(
                                "assay_fetch_failed",
                                extra={
                                    "msg": error_message,
                                    "chunk_ids": context["chunk_ids"],
                                },
                                error=error_message,
                                **log_context,
                            )
                            chunk_failures.add_failure(current, error_message)
                            break

                        delay = compute_backoff_delay(attempt, retry_cfg, jitter=jitter)
                        logger.warning(
                            "assay_fetch_retry",
                            extra={
                                "msg": error_message,
                                "chunk_ids": context["chunk_ids"],
                            },
                            delay=delay,
                            **log_context,
                        )
                        if delay > 0:
                            sleep(delay)
                    else:
                        frames.append(frame)
                        break

            if frames:
                return pd.concat(frames, ignore_index=True)

            return pd.DataFrame(columns=ASSAY_COLUMNS)

        fetch_config = ChunkedFetchConfig(
            ids=ids_source,
            chunk_size=cfg.assay.batch_size,
            workers=1,
        )

        writer_config = CsvWriterConfig(
            writer=write_csv_chunks_deterministic,
            kwargs={
                "chunksize": cfg.io.csv_chunksize,
                "sort_chunksize": cfg.io.csv_chunksize,
                "sep": cfg.io.csv_sep,
                "encoding": cfg.io.csv_encoding,
                "cfg": cfg,
            },
        )

        fetcher, writer = prepare_chunked_pipeline(
            fetch_config=fetch_config,
            fetch_chunk=fetch_chunk,
            csv_writer=writer_config,
        )

        pipeline_stats: dict[str, object] | None = None

        def _capture_stats(stats: Mapping[str, object]) -> None:
            nonlocal pipeline_stats
            pipeline_stats = dict(stats)

        execution = None
        try:
            definition = PipelineDefinition(
                schema=AssaysSchema,
                schema_name="AssaysSchema",
                validators=validators,
                metadata_hooks=metadata_hooks,
                writer=writer,
                command=" ".join(sys.argv),
                config_snapshot=_serialize_paths(cfg.to_dict()),
                inputs={"input_csv": str(args.input_csv)},
                key_columns=["assay_chembl_id"],
                table_quality=table_quality,
                stats_extra=chunk_failures.stats,
                stats_callback=_capture_stats,
                dictionary_resources=dictionary_resources,
            )
            execution = run_pipeline(
                definition=definition,
                fetcher=fetcher,
                output_path=output_path,
                failure_path=failure_path,
                cfg=cfg,
                logger=logger,
                emit_legacy_artifacts=emit_legacy,
            )
        finally:
            if emit_legacy:
                chunk_failures.save(fetch_failure_path, cfg=cfg)
            else:
                Path(fetch_failure_path).unlink(missing_ok=True)
                Path(f"{fetch_failure_path}.meta.yaml").unlink(missing_ok=True)

    if execution is None:
        exit_code_int = 1
        dataset_path: Path | None = None
    else:
        exit_code_attr = getattr(execution, "exit_code", None)
        exit_code_int = int(exit_code_attr if exit_code_attr is not None else execution)
        dataset_path = getattr(execution, "dataset_path", None)

    if exit_code_int == 0 and dataset_path is not None:
        standard_artifacts = _persist_standard_outputs(dataset_path)
        standard_dataset = standard_artifacts.dataset
        if not emit_legacy and standard_dataset != dataset_path:
            Path(dataset_path).unlink(missing_ok=True)
        dataset_path = standard_dataset

    if dataset_path is not None:
        output_path = Path(dataset_path)
        args.final_out = output_path
        args.output_csv = output_path

    exit_status = exit_code_int

    dropped_columns_report = [
        column for column in ASSAY_OUTPUT_DROP_COLUMNS if column in dropped_columns_seen
    ]
    if dropped_columns_report:
        logger.info(
            "Dropped columns from output.assay_*: %s",
            ", ".join(dropped_columns_report),
        )
    else:
        logger.info("Dropped columns from output.assay_*: <none>")

    if limit is not None:
        logger.info("process_limit", limit=processed_ids)
    else:
        logger.info("processed_count", count=processed_ids)

    if pipeline_stats is not None:
        logger.info(
            "records_dropped",
            rows_total=int(pipeline_stats.get("rows_total", processed_ids)),
            rows_kept=int(pipeline_stats.get("rows_kept", 0)),
            rows_dropped=int(pipeline_stats.get("rows_dropped", 0)),
        )

    if exit_status == 0:
        if not postprocess_enabled:
            logger.info(
                "postprocess_skipped",
                output=str(output_path),
                processed=processed_ids,
                message="Postprocessing skipped (--postprocess flag disabled).",
            )
            return exit_status

        assert postprocess_runtime is not None
        assert postprocess_runner is not None
        assert postprocess_generate_metrics is not None
        assert postprocess_pipeline_config is not None
        assert postprocess_csv_cfg is not None
        assert SchemaValidationError is not None
        assert StepError is not None

        output_postprocessed = output_path.with_name(
            "output_postprocessed.assays.csv"
        )

        try:
            postprocess_result = run_postprocessing_pipeline(
                "assays",
                output_path,
                output_postprocessed,
                postprocess_runtime,
            )
        except io.CsvReadError as exc:
            logger.error(
                "assays_postprocess_input_read_failed",
                input=str(output_path),
                error=str(exc.original_error),
            )
            logger.error(
                "assay_pipeline_failed",
                output_postprocessed=str(output_postprocessed),
                processed=processed_ids,
                exit_code=1,
            )
            return 1
        except (SchemaValidationError, StepError) as exc:  # type: ignore[misc]
            logger.exception("assays_postprocess_failed", exc=exc)
            logger.error(
                "assay_pipeline_failed",
                output_postprocessed=str(output_postprocessed),
                processed=processed_ids,
                exit_code=1,
            )
            return 1
        except FileNotFoundError:
            logger.error(
                "assays_postprocess_input_missing",
                input=str(output_path),
            )
            logger.error(
                "assay_pipeline_failed",
                output_postprocessed=str(output_postprocessed),
                processed=processed_ids,
                exit_code=1,
            )
            return 1
        except Exception as exc:  # pragma: no cover - defensive guard
            logger.exception("assays_postprocess_unexpected_error", exc=exc)
            logger.error(
                "assay_pipeline_failed",
                output_postprocessed=str(output_postprocessed),
                processed=processed_ids,
                exit_code=1,
            )
            return 1

        metrics = postprocess_result.metrics
        if metrics is not None:
            summary = metrics.summary()
            summary["output_postprocessed"] = str(postprocess_result.output_path)
            logger.info("assays_postprocess_summary", **summary)

        extras = {
            "input": str(output_path),
            "output_postprocessed": str(postprocess_result.output_path),
            "processed": processed_ids,
        }
        pipeline_version = (
            metrics.pipeline_version
            if metrics and metrics.pipeline_version is not None
            else postprocess_runtime.pipeline_version
            or postprocess_pipeline_config.pipeline_version
            or get_pipeline_version()
        )

        pipeline_metrics, report_path = postprocess_generate_metrics(
            "assays",
            postprocess_result.output_path,
            postprocess_csv_cfg,
            postprocess_runner,
            pipeline_version=pipeline_version,
            extras=extras,
            logger=logger,
            pipeline_metrics=metrics,
        )

        payload: dict[str, Any] = {
            "output_postprocessed": str(postprocess_result.output_path),
            "processed": processed_ids,
            "pipeline_version": pipeline_version,
        }
        if pipeline_metrics is not None:
            summary = pipeline_metrics.summary()
            if summary.get("rows") is not None:
                payload["postprocess_rows"] = summary["rows"]
            if summary.get("columns") is not None:
                payload["postprocess_columns"] = summary["columns"]
            if summary.get("duration_s") is not None:
                payload["postprocess_duration_s"] = summary["duration_s"]
            if summary.get("steps") is not None:
                payload["postprocess_steps"] = summary["steps"]
            if pipeline_metrics.validation is not None:
                payload["postprocess_schema"] = pipeline_metrics.validation.schema
        if report_path is not None:
            payload["postprocess_report"] = str(report_path)
        logger.info("assay_pipeline_done", **payload)
    else:
        logger.error(
            "assay_pipeline_failed",
            output_postprocessed=str(output_path),
            processed=processed_ids,
            exit_code=exit_status,
        )

    return exit_status


def _update_assay_config_from_options(
    cfg: Config, options: AssayPipelineOptions
) -> None:
    """Apply programmatic overrides from ``options`` to ``cfg``."""

    pipelines = cfg.sources.chembl.pipelines
    section = pipelines.assay
    updates: dict[str, object] = {"offset": options.offset}
    if options.limit is not None:
        updates["limit"] = options.limit
    if options.timeout is not None:
        updates["timeout"] = options.timeout
    if options.batch_size is not None:
        updates["batch_size"] = options.batch_size
    pipelines.assay = section.model_copy(update=updates)


def run_assay_service(
    config: Config, options: AssayPipelineOptions
) -> PipelineRunResult:
    """Execute the assay pipeline using typed ``options``."""

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
    _update_assay_config_from_options(cfg, options)

    args = argparse.Namespace(
        input_csv=Path(options.input_csv),
        final_out=output_path,
        output_csv=output_path,
        limit=options.limit,
        offset=options.offset,
        timeout=options.timeout or cfg.assay.timeout,
        batch_size=options.batch_size or cfg.assay.batch_size,
        skip_existing=options.skip_existing,
        force=options.force,
    )

    exit_code = run(cfg, args)
    reason = None if exit_code == 0 else "pipeline_failed"
    written = None if exit_code != 0 else True
    return PipelineRunResult(
        exit_code=exit_code,
        output_path=output_path,
        executed=True,
        reason=reason,
        written=written,
    )


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the assay pipeline handling ``--skip-existing`` semantics."""

    args.emit_legacy_artifacts = bool(getattr(args, "emit_legacy_artifacts", False))

    final_out_attr = getattr(args, "final_out", None)
    if final_out_attr in (None, argparse.SUPPRESS):
        legacy_output = getattr(args, "output_csv", None)
        if legacy_output not in (None, argparse.SUPPRESS):
            output_path = Path(legacy_output)
            if not isinstance(legacy_output, Path):
                args.final_out = output_path
            args.output_csv = output_path
        else:
            output_path = Path(
                io.default_output_path(
                    args.input_csv,
                    cfg.io,
                    date=getattr(args, "date", None),
                )
            )
            args.final_out = output_path
            args.output_csv = output_path
    else:
        output_path = Path(final_out_attr)
        if not isinstance(final_out_attr, Path):
            args.final_out = output_path
        args.output_csv = output_path
    if args.skip_existing and output_path.exists() and not args.force:
        logger.info("pipeline_skip_existing", output=str(output_path))
        return 0
    result = run_chembl(cfg, args)
    return 0 if result is None else int(result)


def _build_parser_impl() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser.

    Returns
    -------
    tuple[argparse.ArgumentParser, LoggerConfig]
        A tuple containing the fully configured parser and default logging
        configuration for the pipeline run.
    """
    parser, log_cfg = base_parser(
        "ChEMBL assay data utilities",
        column="assay_chembl_id",
        chunk_size=10,
        size_option="--batch-size",
        size_dest="batch_size",
        emit_legacy_help=(
            "Write the legacy CSV, metadata and diagnostics alongside the "
            "standard output bundle (default: disabled)."
        ),
    )
    parser.set_defaults(input_csv=Path(DEFAULT_INPUT_NAME))
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help=("Maximum number of identifiers to process; use 0 to skip processing"),
    )
    parser.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    parser.add_argument(
        "--postprocess",
        dest="postprocess",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Enable assay postprocessing after the main pipeline",
    )
    legacy_option = parser._option_string_actions.get("--emit-legacy-artifacts")
    if legacy_option is not None:
        legacy_option.help = (
            "Write the legacy CSV, metadata and diagnostics alongside the standard "
            "output bundle (default: disabled)."
        )
    parser.set_defaults(func=run_chembl)
    return parser, log_cfg


class AssayPipelineCLI(PipelineCLIBase):
    """CLI adapter delegating to the existing assay pipeline helpers."""

    def build_parser(self) -> tuple[argparse.ArgumentParser, LoggerConfig]:
        return _build_parser_impl()

    def prepare_arguments(
        self,
        parser: argparse.ArgumentParser,
        args: argparse.Namespace,
        argv: Sequence[str] | None,
    ) -> argparse.Namespace:
        del parser, argv
        cli.prepare_io_paths(
            args,
            input_default=DEFAULT_INPUT_NAME,
            output_stem=DEFAULT_OUTPUT_STEM,
        )
        return args

    def handle_pre_run(
        self, parser: argparse.ArgumentParser, args: argparse.Namespace
    ) -> int | None:
        if args.limit == 0:
            logger.info("pipeline_skip_limit", limit=args.limit)
            return 0
        if args.limit is not None and args.limit < 0:
            parser.error("--limit must be zero or a positive integer")
        if args.offset < 0:
            parser.error("--offset must be zero or a positive integer")
        return None

    def get_program_name(self) -> str:
        return Path(__file__).with_suffix("").name

    def get_logger(self) -> Logger:
        return logger

    def get_config_mapping(self) -> Mapping[str, str]:
        return {
            "timeout": "assay.timeout",
            "column": "assay.column",
            "batch_size": "assay.batch_size",
            "limit": "assay.limit",
        }

    def run_pipeline(self, cfg: Config, args: argparse.Namespace) -> int:
        return run(cfg, args)


_CLI = AssayPipelineCLI()


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Expose parser creation for legacy tests importing the function."""

    return _CLI.build_parser()


def main(argv: Sequence[str] | None = None) -> int:
    """Delegate to :class:`AssayPipelineCLI` for backwards compatibility."""

    return _CLI.main(argv)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

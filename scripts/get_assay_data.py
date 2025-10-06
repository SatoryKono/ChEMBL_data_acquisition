"""Command line interface for retrieving ChEMBL assay data.

The script provides a reusable :func:`main` entry point together with helper
functions that unit tests import directly. Functions return explicit exit codes
instead of terminating the interpreter to make orchestration easier.
"""

from __future__ import annotations

import sys
from pathlib import Path
from time import sleep

import argparse
from collections import deque
from collections.abc import Iterable, Iterator, Mapping, Sequence
from functools import partial
from itertools import islice

import pandas as pd
import requests

_PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from library.utils.bootstrap import ensure_project_root


if __package__ in {None, ""}:
    ensure_project_root()

from library.integration import chembl_library as cl
from library.pipelines.assay import postprocessing as ap
from library.postprocessing import enrich_assay_metadata
from library import cli
from library import io
from library.common.csv_utils import write_csv_chunks_deterministic
from library.pipelines.assay.chembl_assay import ASSAY_COLUMNS, MAX_ASSAY_CHUNK_SIZE
from library.clients import ChemblClient
from library.common.rate_limiter import get_global_limiter
from library.cli import (
    LoggerConfig,
    ConfigMetadata,
    configure_logger,
)
from library.cli import build_parser as base_parser
from library.cli_utils import run_cli_command, run_pipeline
from library.cli.logging import setup_cli_logging
from library.cli.metadata import prepare_option
from library.config.loader import _serialize_paths
from library.config.models import Config
from library.common.log import logger
from library.pipelines.common import add_pipeline_metadata
from library.table_quality import analyze_table_quality
from library.validation import validate_assays
from library.schemas import AssaysSchema, normalize_assays
from library.pipelines.common import (
    ChunkedFetchConfig,
    CsvWriterConfig,
    prepare_chunked_pipeline,
)
from library.common.fetch_retry import ChunkFailureTracker, compute_backoff_delay

configure_logger = cli.configure_logger

__all__ = ["ap", "configure_logger", "main", "run", "run_chembl"]


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

    allowed_cols = [column for column in df.columns if column not in ASSAY_OUTPUT_DROP_COLUMNS]
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
            setattr(args, "output_csv", output_path)
        else:
            output_path = Path(io.default_output_path(args.input_csv, cfg.io))
            args.final_out = output_path
            setattr(args, "output_csv", output_path)
    else:
        output_path = Path(final_out_attr)
        if not isinstance(final_out_attr, Path):
            args.final_out = output_path
        setattr(args, "output_csv", output_path)
    metadata_obj = getattr(args, "_config_metadata", None)
    if not isinstance(metadata_obj, ConfigMetadata):
        metadata_obj = None
    output_source = "cli" if getattr(args, "final_out", None) else "derived"
    logger.info(
        "assay_pipeline_start",
        input=prepare_option(metadata_obj, value=str(args.input_csv), default_source="cli"),
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
    fetch_failure_path = output_path.with_name(
        f"{output_path.stem}_fetch_failures.csv"
    )

    def _enrich_with_dictionaries(frame: pd.DataFrame) -> pd.DataFrame:
        return enrich_assay_metadata(frame, dictionary_dir=cfg.resources.dictionary_dir)

    dropped_columns_seen: set[str] = set()

    def _drop_output_columns(frame: pd.DataFrame) -> pd.DataFrame:
        removed = [column for column in ASSAY_OUTPUT_DROP_COLUMNS if column in frame.columns]
        if removed:
            dropped_columns_seen.update(removed)
        return remove_assay_output_columns(frame)

    metadata_hooks = [
        _enrich_with_dictionaries,
        ap.postprocess_assays,
        normalize_assays,
        add_pipeline_metadata,
        _drop_output_columns,
        _drop_assay_output_columns,
    ]

    validators = [partial(validate_assays, return_result=True)]

    doc_quality_cfg = cfg.system.doc_quality
    if doc_quality_cfg.enable:
        table_quality = partial(
            analyze_table_quality,
            table_name=str(output_path.with_suffix("")),
            destination_dir=output_path.parent,
            sample_rows=doc_quality_cfg.sample_rows,
            include_columns=doc_quality_cfg.include_columns,
            exclude_columns=doc_quality_cfg.exclude_columns,
        )
    else:
        def table_quality(_: Path) -> None:
            return None

    rate_cfg = cfg.rate
    global_limiter = None
    if (rate_cfg.global_rps or 0) > 0:
        global_limiter = get_global_limiter(rate_cfg.global_rps, rate_cfg.global_burst)

    with ChemblClient(
        cfg.api,
        cfg.retry,
        cfg.chembl,
        global_limiter=global_limiter,
    ) as client:

        retry_cfg = cfg.retry
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
                        log_context = {k: v for k, v in context.items() if k != "chunk_ids"}

                        if attempt >= attempts:
                            if len(current) > 1:
                                split_index = max(1, len(current) // 2)
                                left = current[:split_index]
                                right = current[split_index:]
                                logger.warning(
                                    "assay_fetch_split",
                                    extra={"msg": error_message, "chunk_ids": context["chunk_ids"]},
                                    **log_context,
                                )
                                if right:
                                    pending.appendleft(right)
                                if left:
                                    pending.appendleft(left)
                                break

                            logger.error(
                                "assay_fetch_failed",
                                extra={"msg": error_message, "chunk_ids": context["chunk_ids"]},
                                error=error_message,
                                **log_context,
                            )
                            chunk_failures.add_failure(current, error_message)
                            break

                        delay = compute_backoff_delay(attempt, retry_cfg)
                        logger.warning(
                            "assay_fetch_retry",
                            extra={"msg": error_message, "chunk_ids": context["chunk_ids"]},
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

        try:
            exit_code = run_pipeline(
                fetcher=fetcher,
                schema=AssaysSchema,
                schema_name="AssaysSchema",
                validators=validators,
                metadata_hooks=metadata_hooks,
                writer=writer,
                output_path=output_path,
                failure_path=failure_path,
                command=" ".join(sys.argv),
                config_snapshot=_serialize_paths(cfg.to_dict()),
                inputs={"input_csv": str(args.input_csv)},
                key_columns=["assay_chembl_id"],
                table_quality=table_quality,
                cfg=cfg,
                stats_extra=chunk_failures.stats,
                logger=logger,
                stats_callback=_capture_stats,
            )
        finally:
            chunk_failures.save(fetch_failure_path, cfg=cfg)

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

    if exit_code == 0:
        logger.info(
            "assay_pipeline_done",
            output=str(output_path),
            processed=processed_ids,
        )
    else:
        logger.error(
            "assay_pipeline_failed",
            output=str(output_path),
            processed=processed_ids,
            exit_code=exit_code,
        )

    return exit_code


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the assay pipeline handling ``--skip-existing`` semantics."""

    final_out_attr = getattr(args, "final_out", None)
    if final_out_attr in (None, argparse.SUPPRESS):
        legacy_output = getattr(args, "output_csv", None)
        if legacy_output not in (None, argparse.SUPPRESS):
            output_path = Path(legacy_output)
            if not isinstance(legacy_output, Path):
                args.final_out = output_path
            setattr(args, "output_csv", output_path)
        else:
            output_path = Path(io.default_output_path(args.input_csv, cfg.io))
            args.final_out = output_path
            setattr(args, "output_csv", output_path)
    else:
        output_path = Path(final_out_attr)
        if not isinstance(final_out_attr, Path):
            args.final_out = output_path
        setattr(args, "output_csv", output_path)
    if args.skip_existing and output_path.exists() and not args.force:
        logger.info("pipeline_skip_existing", output=str(output_path))
        return 0
    return run_chembl(cfg, args)


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
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
        help=(
            "Maximum number of identifiers to process; use 0 to skip processing"
        ),
    )
    parser.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    parser.set_defaults(func=run_chembl)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Execute the assay pipeline with optional argument overrides.

    Parameters
    ----------
    argv : Sequence[str] | None, optional
        Command-line arguments to parse. When ``None`` the process arguments
        from :data:`sys.argv` are used.

    Returns
    -------
    int
        ``0`` on success, non-zero otherwise. Failures are logged with context
        describing the failing section.

    Raises
    ------
    SystemExit
        Raised when invalid command-line options are supplied.
    """
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    cli.prepare_io_paths(
        args,
        input_default=DEFAULT_INPUT_NAME,
        output_stem=DEFAULT_OUTPUT_STEM,
    )
    if args.limit == 0:
        logger.info("pipeline_skip_limit", limit=args.limit)
        return 0
    if args.limit is not None and args.limit < 0:
        parser.error("--limit must be zero or a positive integer")
    if args.offset < 0:
        parser.error("--offset must be zero or a positive integer")
    with setup_cli_logging(
        Path(__file__).with_suffix("").name, log_cfg, getattr(args, "date", None)
    ) as logging_ctx:
        exit_code = run_cli_command(
            args=args,
            parser=parser,
            log_cfg=logging_ctx.log_cfg,
            mapping={
                "timeout": "assay.timeout",
                "column": "assay.column",
                "batch_size": "assay.batch_size",
                "limit": "assay.limit",
            },
            run=run,
            logger=logger,
        )
    configure_logger(log_cfg)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

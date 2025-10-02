"""Command line interface for retrieving ChEMBL assay data.

The script provides a reusable :func:`main` entry point together with helper
functions that unit tests import directly. Functions return explicit exit codes
instead of terminating the interpreter to make orchestration easier.
"""

from __future__ import annotations

import sys
from pathlib import Path

import argparse
from collections.abc import Iterable, Iterator, Sequence
from functools import partial
from itertools import islice

import pandas as pd
import requests

from library.utils.bootstrap import ensure_project_root


if __package__ in {None, ""}:
    ensure_project_root()

from library import assay_postprocessing as ap
from library import chembl_library as cl
from library import cli
from library import io
from library.csv_utils import write_csv_chunks_deterministic
from library.clients import ChemblClient
from library.rate_limiter import get_global_limiter
from library.cli import (
    LoggerConfig,
)
from library.cli import build_parser as base_parser
from library.cli_utils import PipelineError, run_cli_command, run_pipeline
from library.config import Config, _serialize_paths
from library.log import logger
from library.pipeline_metadata import add_pipeline_metadata
from library.table_quality import analyze_table_quality
from library.validation import validate_assays
from schemas import AssaysSchema, normalize_assays

__all__ = ["ap", "main", "run", "run_chembl"]


DEFAULT_INPUT_NAME = "assay.csv"
DEFAULT_OUTPUT_STEM = "assays"


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute assay retrieval from the ChEMBL API.

    The output CSV arranges columns so that fields defined in
    :class:`~schemas.assays.AssaysSchema` appear first. Any additional columns
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

    output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
    failure_path = Path(output).with_name(f"{Path(output).stem}_failure_cases.csv")

    def fetcher() -> Iterable[pd.DataFrame]:
    global_limiter = get_global_limiter(cfg.rate.global_rps, cfg.rate.global_burst)

    with ChemblClient(
        cfg.api, cfg.retry, cfg.chembl, global_limiter=global_limiter
    ) as client:
            chunk_iter = cl._chunked(ids_source, cfg.assay.batch_size)
            for chunk_ids in chunk_iter:
                try:
                    df = cl.get_assays(
                        chunk_ids,
                        cfg=cfg.api,
                        client=client,
                        chunk_size=cfg.assay.batch_size,
                        timeout=cfg.assay.timeout,
                    )
                except (requests.RequestException, ValueError) as exc:
                    logger.error(
                        "assay_fetch_failed",
                        extra={"msg": str(exc)},
                        error=str(exc),
                        batch_size=cfg.assay.batch_size,
                        timeout=cfg.assay.timeout,
                    )
                    raise PipelineError(str(exc)) from exc
                yield df

    metadata_hooks = [
        ap.postprocess_assays,
        normalize_assays,
        add_pipeline_metadata,
    ]

    validators = [partial(validate_assays, return_result=True)]

    def writer(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        col_order: Sequence[str],
        key_cols: Sequence[str],
    ) -> Path:
        sort_columns = list(key_cols) or sorted(col_order)
        return write_csv_chunks_deterministic(
            chunks,
            destination,
            key_cols=sort_columns,
            col_order=col_order,
            chunksize=cfg.io.csv_chunksize,
            sort_chunksize=cfg.io.csv_chunksize,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            cfg=cfg,
        )

    table_quality = partial(
        analyze_table_quality,
        table_name=str(Path(output).with_suffix("")),
    )

    exit_code = run_pipeline(
        fetcher=fetcher,
        schema=AssaysSchema,
        schema_name="AssaysSchema",
        validators=validators,
        metadata_hooks=metadata_hooks,
        writer=writer,
        output_path=output,
        failure_path=failure_path,
        command=" ".join(sys.argv),
        config_snapshot=_serialize_paths(cfg.to_dict()),
        inputs={"input_csv": str(args.input_csv)},
        key_columns=["assay_chembl_id"],
        table_quality=table_quality,
        cfg=cfg,
        logger=logger,
    )

    if limit is not None:
        logger.info("process_limit", limit=processed_ids)
    else:
        logger.info("processed_count", count=processed_ids)

    return exit_code


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the assay pipeline handling ``--skip-existing`` semantics."""

    output_path = Path(
        args.output_csv or io.default_output_path(args.input_csv, cfg.io)
    )
    args.output_csv = output_path
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
    return run_cli_command(
        args=args,
        parser=parser,
        log_cfg=log_cfg,
        mapping={
            "timeout": "assay.timeout",
            "column": "assay.column",
            "batch_size": "assay.batch_size",
            "limit": "assay.limit",
        },
        run=run,
        logger=logger,
    )


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

"""Command line interface for retrieving ChEMBL activity data."""

from __future__ import annotations

import argparse
import sys
from typing import Sequence

import requests
from library.config import Config, ensure_dirs, print_config, _serialize_paths
from library.chembl_client import init_session

from library import chembl_library as cl
from library import io
from library.metadata import Stats, file_sha256, write_meta_yaml
from library.cli import (
    apply_config_overrides,
    build_parser as base_parser,
    configure_logger,
    LoggerConfig,
)
from library.log import logger
from library.sidecar import SidecarErrors
from library.table_quality import analyze_table_quality
from pandera.errors import SchemaErrors
from schemas import ActivitiesSchema, normalize_activities

from library import write_csv_deterministic

ORIG_WRITE_CSV = io.write_csv


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute activity retrieval from the ChEMBL API.

    Parameters
    ----------
    cfg : Config
        Application configuration.
    args : argparse.Namespace
        Parsed command-line arguments. ``args.limit`` constrains the number of
        identifiers processed, while ``args.dry_run`` skips network calls and
        file generation.

    Returns
    -------
    int
        Zero on success, non-zero on failure.

    """
    if args.limit is not None and args.limit < 0:
        logger.error("--limit must be non-negative")
        return 1

    if args.dry_run:
        expected = args.limit if args.limit is not None else 0
        logger.info("dry run selected; would process at most %d identifiers", expected)
        return 0

    # Configure HTTP session with the supplied User-Agent and retry policy
    init_session(cfg.api, cfg.retry)

    try:
        ids = io.read_ids(
            args.input_csv,
            column=args.column,
            cfg=cfg.io,
            sep=args.sep,
            encoding=args.encoding,
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error("%s", exc)
        return 1

    if args.limit is not None:
        ids = ids[: args.limit]
        logger.info("processing at most %d identifiers", len(ids))

    try:
        df = cl.get_activities(
            ids, cfg=cfg.api, chunk_size=args.chunk_size, timeout=args.timeout
        )
    except (requests.RequestException, ValueError) as exc:
        logger.error("failed to retrieve activities: %s", exc)
        return 1
    output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
    df = normalize_activities(df)
    rows_total = len(df)
    exit_code = 0
    required_cols = set(ActivitiesSchema.columns.keys())
    if required_cols.issubset(df.columns):
        try:
            df = ActivitiesSchema.validate(df, lazy=True)
        except SchemaErrors as exc:
            failure_path = output.with_name(f"{output.stem}_failure_cases.csv")
            errors = SidecarErrors()
            for row in exc.failure_cases.to_dict("records"):
                errors.add_error(row)
            errors.save(failure_path)
            logger.error(
                "validation failed; wrote %d failure cases to %s",
                len(exc.failure_cases),
                failure_path,
            )
            df = getattr(exc, "validated_data", df)
            exit_code = 1
    else:
        missing = required_cols.difference(df.columns)
        logger.warning("Skipping validation due to missing columns: %s", missing)
    rows_kept = len(df)
    rows_dropped = rows_total - rows_kept
    try:
        key_cols = [c for c in ["activity_id"] if c in df.columns]
        csv_path = write_csv_deterministic(
            df,
            output,
            col_order=[
                "activity_id",
                "testitem_id",
                "target_id",
                "standard_type",
                "standard_value",
                "pA_value",
            ],
            key_cols=key_cols or None,
        )
        if io.write_csv is not ORIG_WRITE_CSV:
            io.write_csv(df, output, cfg=cfg, sep=args.sep, encoding=args.encoding)
        logger.info("Wrote %d rows to %s", rows_kept, csv_path)
    except OSError as exc:
        logger.error("failed to write output CSV: %s", exc)
        return 1

    stats: Stats = {
        "rows_total": rows_total,
        "rows_kept": rows_kept,
        "rows_dropped": rows_dropped,
        "output_sha256": file_sha256(csv_path),
    }
    write_meta_yaml(
        csv_path=csv_path,
        command=" ".join(sys.argv),
        config_subset=_serialize_paths(cfg.to_dict()),
        inputs={"input_csv": str(args.input_csv)},
        stats=stats,
        schema="ActivitiesSchema",
    )

    try:
        analyze_table_quality(df, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return exit_code


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser."""
    parser, log_cfg = base_parser(
        "ChEMBL activity data utilities", column="activity_id", chunk_size=5
    )
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
        help="Maximum number of identifiers to process",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Read input and exit without contacting the API or writing files",
    )
    parser.set_defaults(func=run_chembl)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults."""
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline start run_id=%s", log_cfg.run_id, extra={"event": "start"})
    try:
        cfg: Config = apply_config_overrides(
            args,
            parser,
            args.config,
            mapping={"timeout": "api.timeout_read"},
        )
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
            logger.info(
                "pipeline done run_id=%s", log_cfg.run_id, extra={"event": "done"}
            )
            return 0
        ensure_dirs(cfg)
        logger = configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
    except (ValueError, TypeError) as exc:
        logger.error("%s", exc)
        logger.info("pipeline fail run_id=%s", log_cfg.run_id, extra={"event": "fail"})
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error("failed to set up directories: %s", exc)
        logger.info("pipeline fail run_id=%s", log_cfg.run_id, extra={"event": "fail"})
        return 1
    exit_code = args.func(cfg, args)
    if exit_code == 0:
        logger.info("pipeline done run_id=%s", log_cfg.run_id, extra={"event": "done"})
    else:
        logger.info("pipeline fail run_id=%s", log_cfg.run_id, extra={"event": "fail"})
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

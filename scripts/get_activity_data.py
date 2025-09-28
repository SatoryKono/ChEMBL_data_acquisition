"""Command line interface for retrieving ChEMBL activity data."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    project_root = Path(__file__).resolve().parents[1]
    project_root_str = str(project_root)
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)
    import bootstrap  # type: ignore[import-not-found]
else:
    from . import bootstrap  # type: ignore[import-not-found]

bootstrap.ensure_project_root()

from library import cli, io
from library.cli import LoggerConfig, configure_logger
from library.cli import build_parser as base_parser
from library.log import logger
from library.processing.activity import (
    ActivityFetchError,
    ActivityReadError,
    fetch_activity_frame,
    write_activity_outputs,
)
from library.utils.config import Config, ConfigError, ensure_dirs, print_config


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute activity retrieval from the ChEMBL API."""

    limit = cfg.activity.limit
    if limit is not None and limit < 0:
        logger.error("invalid_limit", section="activity.limit", limit=limit)
        return 1

    if cfg.activity.dry_run:
        expected = limit if limit is not None else 0
        logger.info("dry_run", limit=expected)
        return 0

    try:
        frame = fetch_activity_frame(cfg, args.input_csv, limit=limit)
    except ActivityReadError as exc:
        logger.error(
            "read_fail",
            error=str(exc.error),
            path=str(exc.path),
        )
        return 1
    except ActivityFetchError as exc:
        logger.error(
            "activity_fetch_failed",
            extra={"msg": str(exc.error)},
            error=str(exc.error),
            batch_size=cfg.activity.batch_size,
            timeout=cfg.activity.timeout,
        )
        return 1

    output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
    try:
        result = write_activity_outputs(
            frame,
            output,
            cfg,
            command=" ".join(sys.argv),
            input_csv=args.input_csv,
        )
    except OSError as exc:
        logger.error(
            "write_fail",
            error=str(exc),
            path=str(output),
        )
        return 1
    except ValueError as exc:
        logger.error(
            "quality_report_failed",
            error=str(exc),
            path=str(output),
        )
        return 1

    return result.exit_code


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser."""

    parser, log_cfg = base_parser(
        "ChEMBL activity data utilities",
        column="activity_id",
        chunk_size=5,
        size_option="--batch-size",
        size_dest="batch_size",
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
    if args.limit is not None and args.limit <= 0:
        parser.error("--limit must be a positive integer")
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline_start", run_id=log_cfg.run_id)
    try:
        cfg: Config = cli.apply_config_overrides(
            args,
            parser,
            args.config,
            mapping={
                "timeout": "activity.timeout",
                "column": "activity.column",
                "batch_size": "activity.batch_size",
                "limit": "activity.limit",
                "dry_run": "activity.dry_run",
            },
        )
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        logger = configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
    except (ConfigError, ValueError, TypeError) as exc:
        logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error(
            "directory_setup_failed",
            error=str(exc),
            output=str(args.output_csv),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    exit_code: int = args.func(cfg, args)
    if exit_code == 0:
        logger.info("pipeline_done", run_id=log_cfg.run_id)
    else:
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

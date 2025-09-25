"""CLI entry point for table quality analysis."""

from __future__ import annotations

import sys

# ruff: noqa: E402
from pathlib import Path

if __package__ is None:  # running as a script
    sys.path.append(str(Path(__file__).resolve().parents[1]))

import argparse
import os
from collections.abc import Sequence

import pandas as pd

from library.cli import (
    LoggerConfig,
    apply_config_overrides,
    configure_logger,
)
from library.cli import (
    build_parser as base_parser,
)
from library.config import Config, ensure_dirs, print_config
from library.log import logger
from library.table_quality import analyze_table_quality


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the profiling workflow.

    Parameters
    ----------
    cfg : Config
        Application configuration.
    args:
        Parsed command-line arguments.

    Returns
    -------
    int
        ``0`` on success, non-zero on failure.

    """
    try:
        try:
            df = pd.read_csv(args.input_csv, sep=args.sep, encoding=args.encoding)
        except (FileNotFoundError, pd.errors.ParserError, UnicodeError) as exc:
            logger.error(
                "input_read_failed",
                error=str(exc),
                path=str(args.input_csv),
                encoding=args.encoding,
                sep=args.sep,
            )
            return 1

        original_cwd = Path.cwd()
        try:
            args.output_csv.mkdir(parents=True, exist_ok=True)
            os.chdir(args.output_csv)
            analyze_table_quality(df, table_name=args.table_name)
        finally:
            os.chdir(original_cwd)
        return 0
    except Exception as exc:  # pragma: no cover - defensive
        logger.exception("run_fail", exc=exc)
        return 1


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser."""
    parser, log_cfg = base_parser("Table quality analysis", column="chembl_id")
    parser.add_argument(
        "--table-name",
        required=True,
        help="Base name used for output report files",
    )
    parser.set_defaults(func=run, output_csv=Path("."), encoding="utf-8-sig")
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults."""
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline_start", run_id=log_cfg.run_id)
    try:
        cfg: Config = apply_config_overrides(args, parser, args.config)
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        logger = configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
    except (ValueError, TypeError) as exc:
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

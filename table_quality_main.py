"""CLI entry point for table quality analysis."""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
from typing import Sequence

import pandas as pd
from library.config import Config, ensure_dirs
from library.cli import add_config_argument, apply_config_overrides, configure_logging

from library.table_quality import analyze_table_quality

logger = logging.getLogger(__name__)


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the profiling workflow.

    Parameters
    ----------
    args:
        Parsed command-line arguments.

    Returns
    -------
    int
        ``0`` on success, non-zero on failure.

    """
    try:
        df = pd.read_csv(args.input_csv, sep=args.sep, encoding=args.encoding)
    except (FileNotFoundError, pd.errors.ParserError, UnicodeError) as exc:
        logger.error("%s", exc)
        return 1

    original_cwd = Path.cwd()
    try:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        os.chdir(args.output_dir)
        analyze_table_quality(df, table_name=args.table_name)
    finally:
        os.chdir(original_cwd)
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""
    parser = argparse.ArgumentParser(description="Table quality analysis")
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    parser.add_argument(
        "--input",
        dest="input_csv",
        type=Path,
        default=Path("input.csv"),
        help="Input CSV file",
    )
    parser.add_argument(
        "--table-name",
        required=True,
        help="Base name used for output report files",
    )
    parser.add_argument(
        "--output",
        dest="output_dir",
        type=Path,
        default=Path("."),
        help="Directory to store generated reports",
    )
    parser.add_argument("--sep", default=",", help="CSV delimiter")
    parser.add_argument("--encoding", default="utf-8-sig", help="File encoding")
    parser.add_argument(
        "--print-config",
        action="store_true",
        help="Print effective configuration and exit",
    )
    parser.set_defaults(func=run)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point."""
    parser = build_parser()
    add_config_argument(parser)
    args = parser.parse_args(argv)
    cfg = apply_config_overrides(args, parser, args.config)
    ensure_dirs(cfg)
    if args.print_config:
        print(cfg.to_yaml())
        return 0
    configure_logging(args.log_level, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
    return args.func(cfg, args)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

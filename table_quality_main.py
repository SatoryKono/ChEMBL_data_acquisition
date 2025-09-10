"""CLI entry point for table quality analysis."""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
from typing import Sequence

import pandas as pd

from library.table_quality import analyze_table_quality

from library.config import DEFAULT_CONFIG load_config


logger = logging.getLogger(__name__)


def run(args: argparse.Namespace) -> int:
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
        "--config",
        type=Path,
        default=Path("config.yaml"),
        help="Path to YAML configuration file",
    )
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

        default=DEFAULT_CONFIG.output.data_dir,

        help="Directory to store generated reports",
    )
    parser.add_argument("--sep", default=",", help="CSV delimiter")
    parser.add_argument("--encoding", default="utf-8-sig", help="File encoding")
    parser.set_defaults(func=run)
    return parser


def configure_logging(level: str) -> None:
    """Configure basic logging."""
    numeric_level = getattr(logging, level.upper(), logging.INFO)
    logging.basicConfig(level=numeric_level)


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    config = load_config(args.config)
    if args.output_dir is None:
        args.output_dir = Path(config.get("output", {}).get("data_dir", "."))
    configure_logging(args.log_level)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

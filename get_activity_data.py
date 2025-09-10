"""Command line interface for retrieving ChEMBL activity data."""

from __future__ import annotations

import argparse
import logging
from typing import Sequence

import requests

from library import chembl_library as cl
from library import io
from library.cli import build_parser as base_parser, configure_logging
from library.table_quality import analyze_table_quality
from library.config import DEFAULT_CONFIG

logger = logging.getLogger(__name__)


def run_chembl(args: argparse.Namespace) -> int:
    """Execute activity retrieval from the ChEMBL API.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    int
        Zero on success, non-zero on failure.

    """
    try:
        ids = io.read_ids(
            args.input_csv,
            column=args.column,
            sep=args.sep,
            encoding=args.encoding,
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error("%s", exc)
        return 1

    try:
        df = cl.get_activities(ids, chunk_size=args.chunk_size, timeout=args.timeout)
    except (requests.RequestException, ValueError) as exc:
        logger.error("failed to retrieve activities: %s", exc)
        return 1
    output = args.output_csv or io.default_output_path(args.input_csv)
    try:
        io.write_csv(df, output, sep=args.sep, encoding=args.encoding)
        logger.info("Wrote %d rows to %s", len(df), output)
    except OSError as exc:
        logger.error("failed to write output CSV: %s", exc)
        return 1
    try:
        analyze_table_quality(df, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""
    parser = base_parser(
        "ChEMBL activity data utilities", column="activity_id", chunk_size=5
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=DEFAULT_CONFIG.timeouts.read,
        help="Timeout in seconds for each HTTP request",
    )
    parser.set_defaults(func=run_chembl)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.log_level)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

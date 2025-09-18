"""Command line interface for ChEMBL activity extraction."""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

from library.activity_extraction import extract_activities
from library.config import Config


def _positive_int(value: str) -> int:
    """Return ``value`` converted to :class:`int` ensuring it is positive."""

    integer = int(value)
    if integer <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return integer


def build_parser() -> argparse.ArgumentParser:
    """Create the argument parser for the activity extraction CLI."""

    parser = argparse.ArgumentParser(
        description="Download ChEMBL activity data into a validated CSV."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("input.csv"),
        help="Path to the CSV file containing activity identifiers.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Destination for the generated CSV file.",
    )
    parser.add_argument(
        "--column",
        default="activity_id",
        help="Column name holding the activity identifiers.",
    )
    parser.add_argument(
        "--chunk-size",
        type=_positive_int,
        default=5,
        help="Maximum number of identifiers requested per API call.",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds applied to API requests.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional cap on the number of identifiers processed.",
    )
    parser.add_argument(
        "--sep",
        default=None,
        help="Field delimiter of the input CSV (defaults to configuration value).",
    )
    parser.add_argument(
        "--encoding",
        default=None,
        help="Encoding of the input CSV (defaults to configuration value).",
    )
    parser.add_argument(
        "--na-marker",
        action="append",
        dest="na_markers",
        help="Additional marker interpreted as a missing value.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging level (e.g. DEBUG, INFO, WARNING).",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Parse arguments and invoke :func:`extract_activities`."""

    parser = build_parser()
    args = parser.parse_args(argv)

    if args.limit is not None and args.limit < 0:
        parser.error("--limit must be non-negative")
    if args.timeout < 0:
        parser.error("--timeout must be non-negative")

    level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(level, int):
        parser.error(f"invalid log level: {args.log_level}")

    logging.basicConfig(
        level=level,
        format="[%(asctime)s] %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    cfg = Config()

    return extract_activities(
        input_csv=args.input,
        output_csv=args.output,
        column=args.column,
        chunk_size=args.chunk_size,
        timeout=args.timeout,
        limit=args.limit,
        sep=args.sep,
        encoding=args.encoding,
        na_markers=args.na_markers,
        cfg=cfg,
        command=" ".join(sys.argv),
    )


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

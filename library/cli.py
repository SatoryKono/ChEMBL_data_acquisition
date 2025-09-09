"""Shared command-line helpers."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path


def build_parser(
    description: str, *, column: str, chunk_size: int = 10
) -> argparse.ArgumentParser:
    """Return a parser with common arguments.

    Parameters
    ----------
    description:
        Text used in the parser description.
    column:
        Default column name for identifier extraction.
    chunk_size:
        Default chunk size for API requests.
    """

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    parser.add_argument(
        "--input",
        dest="input_csv",
        type=Path,
        default=Path("input.csv"),
        help="Input CSV file",
    )
    parser.add_argument(
        "--output",
        dest="output_csv",
        type=Path,
        default=None,
        help="Destination CSV file (default: auto-generate)",
    )
    parser.add_argument(
        "--column", default=column, help="Identifier column in input CSV"
    )
    parser.add_argument("--sep", default=",", help="CSV delimiter")
    parser.add_argument("--encoding", default="utf8", help="File encoding")
    parser.add_argument(
        "--chunk-size", type=int, default=chunk_size, help="Maximum IDs per request"
    )
    return parser


def configure_logging(level: str) -> None:
    """Configure basic logging for CLI entry points."""
    numeric = getattr(logging, level.upper(), logging.INFO)
    logging.basicConfig(level=numeric)

"""Command line interface for retrieving ChEMBL activity data."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path
from typing import Sequence

import pandas as pd

from library import chembl_library as cl

logger = logging.getLogger(__name__)


def read_ids(
    path: str | Path,
    column: str = "activity_id",
    sep: str = ",",
    encoding: str = "utf8",
) -> list[str]:
    """Read ChEMBL activity identifiers from a CSV file.

    Parameters
    ----------
    path : str or Path
        Path to the CSV file.
    column : str, optional
        Name of the column containing activity identifiers. Defaults to
        ``"activity_id"``.
    sep : str, optional
        Field delimiter, by default a comma.
    encoding : str, optional
        File encoding, by default ``"utf8"``.

    Returns
    -------
    list[str]
        Identifier values in the order they appear. Empty strings and
        ``"#N/A"`` markers are discarded.

    Raises
    ------
    ValueError
        If ``column`` is not present in the input file.
    """
    try:
        with Path(path).open("r", encoding=encoding, newline="") as fh:
            reader = csv.DictReader(fh, delimiter=sep)
            if reader.fieldnames is None or column not in reader.fieldnames:
                raise ValueError(f"column '{column}' not found in {path}")
            ids: list[str] = []
            for row in reader:
                value = (row.get(column) or "").strip()
                if value and value != "#N/A":
                    ids.append(value)
            return ids
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"input file not found: {path}") from exc
    except csv.Error as exc:
        raise ValueError(f"malformed CSV in file: {path}: {exc}") from exc


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
        ids = read_ids(
            args.input_csv,
            column=args.column,
            sep=args.sep,
            encoding=args.encoding,
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error("%s", exc)
        return 1

    df = cl.get_activities(ids, chunk_size=args.chunk_size, timeout=args.timeout)
    try:
        df.to_csv(args.output_csv, index=False, sep=args.sep, encoding=args.encoding)
        logger.info("Wrote %d rows to %s", len(df), args.output_csv)
        return 0
    except OSError as exc:
        logger.error("failed to write output CSV: %s", exc)
        return 1


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""
    parser = argparse.ArgumentParser(description="ChEMBL activity data utilities")
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    parser.add_argument(
        "input_csv", type=Path, help="CSV file containing activity identifiers"
    )
    parser.add_argument("output_csv", type=Path, help="Destination CSV file")
    parser.add_argument(
        "--column",
        default="activity_id",
        help="Column name in the input CSV containing identifiers",
    )
    parser.add_argument("--sep", default=",", help="CSV delimiter")
    parser.add_argument("--encoding", default="utf8", help="File encoding")
    parser.add_argument(
        "--chunk-size", type=int, default=5, help="Maximum number of IDs per request"
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    parser.set_defaults(func=run_chembl)
    return parser


def configure_logging(level: str) -> None:
    """Configure basic logging."""
    numeric_level = getattr(logging, level.upper(), logging.INFO)
    logging.basicConfig(level=numeric_level)


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.log_level)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

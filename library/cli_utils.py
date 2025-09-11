"""Argument parser helpers for lightweight command-line interfaces."""

from __future__ import annotations

import argparse

from .cli import add_common_arguments


def build_parser() -> argparse.ArgumentParser:
    """Create an ``argparse`` parser with shared CSV options.

    Returns
    -------
    argparse.ArgumentParser
        Parser pre-configured with common command-line flags used by
        :mod:`csv_utils_main`.
    """
    parser = argparse.ArgumentParser(
        description=(
            "CLI wrapper for :func:`write_csv_deterministic`. This script reads "
            "an input CSV file and re-serialises it deterministically using "
            ":func:`library.csv_utils.write_csv_deterministic`."
        )
    )
    add_common_arguments(parser)
    parser.add_argument("--col-order", nargs="*", help="Preferred column order")
    parser.add_argument("--key-cols", nargs="*", help="Columns used for sorting")
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=1000,
        help="Number of rows read per chunk",
    )
    return parser

"""CLI entry point for the dtype inspector.

The script executes a small sample retrieval for each ``get_*_data`` helper and
logs the resulting pandas dtypes.  Use this during development to spot schema
changes when upstream APIs evolve.

Example
-------
Run with default settings::

    python scripts/dtype_inspector_main.py --log-level INFO
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence

from library.cli import LoggerConfig, configure_logger
from library.dtype_inspector import inspect_dtypes


def build_parser() -> argparse.ArgumentParser:
    """Return the argument parser."""

    parser = argparse.ArgumentParser(description="Inspect dtypes from sample data")
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Program entry point."""

    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logger(LoggerConfig(level=args.log_level))
    inspect_dtypes()
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

"""CLI for combining ChEMBL initialisation Excel sources."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Sequence

from library import input_initialisation_library as lib

logger = logging.getLogger(__name__)


def run(args: argparse.Namespace) -> int:
    """Execute table combination routine.

    Parameters
    ----------
    args:
        Parsed command line arguments.

    Returns
    -------
    int
        Zero on success, non-zero otherwise.
    """
    try:
        if not args.same_doc.exists():
            raise FileNotFoundError(f"{args.same_doc} does not exist")
        if not args.all_doc.exists():
            raise FileNotFoundError(f"{args.all_doc} does not exist")
        args.out_dir.mkdir(parents=True, exist_ok=True)

        logger.info("Loading source workbooks")
        same = lib.load_same_doc(args.same_doc)
        all_ = lib.load_all_doc(args.all_doc)

        logger.info("Combining tables")
        tables = lib.build_combined_tables(same, all_)

        logger.info("Saving output")
        lib.save_tables(tables, args.out_dir, fmt=args.format)
        return 0
    except Exception as exc:
        logger.error("%s", exc)
        return 1


def build_parser() -> argparse.ArgumentParser:
    """Create argument parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    parser.add_argument(
        "--same-doc", type=Path, required=True, help="Path to same document workbook"
    )
    parser.add_argument(
        "--all-doc", type=Path, required=True, help="Path to all document workbook"
    )
    parser.add_argument(
        "--out-dir", type=Path, default=Path("."), help="Output directory"
    )
    parser.add_argument(
        "--format", choices=["csv"], default="csv", help="Output format"
    )
    parser.set_defaults(func=run)
    return parser


def configure_logging(level: str) -> None:
    """Configure logging at ``level``."""
    numeric = getattr(logging, level.upper(), logging.INFO)
    logging.basicConfig(level=numeric)


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.log_level)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

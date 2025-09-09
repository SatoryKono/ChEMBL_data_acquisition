"""CLI for combining ChEMBL initialisation Excel sources.

Exports pair tables without merging:

- ``pairs_same_document.csv`` from sheet ``pairs_same_doc`` in ``--same-doc``

- ``pairs_independent.csv`` and ``pairs_non_independent.csv`` derived from
  ``step5_pairs`` in ``--all-doc`` based on the ``INDEPENDENT`` flag
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Sequence

from library import input_initialisation_library as lib
from library.table_quality import analyze_table_quality

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
        tables = lib.build_combined_tables(
            same, all_, dictionary_dir=args.dictionary_dir
        )

        logger.info("Saving output")
        paths = lib.save_tables(tables, args.out_dir, fmt=args.format)
        # Ensure that files were actually written to disk
        missing = [str(p) for p in paths.values() if not p.exists()]
        if missing:
            raise RuntimeError("failed to write output files: " + ", ".join(missing))

        logger.info("Generating data quality reports")
        report_dir = args.out_dir / "data_validity_report"
        report_dir.mkdir(parents=True, exist_ok=True)
        for entity, path in paths.items():
            logger.info("Profiling %s", entity)
            analyze_table_quality(path, table_name=str(report_dir / path.stem))

        logger.info(
            "Saved %d tables and quality reports to %s", len(paths), args.out_dir
        )
        return 0
    except KeyError as exc:
        logger.error("required table '%s' missing", exc.args[0])
        return 1
    except Exception as exc:  # pragma: no cover - defensive
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
        "--dictionary-dir",
        type=Path,
        default=Path("dictionary"),
        help="Directory with targets_type.csv, citation_fraction.csv and status.csv",
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

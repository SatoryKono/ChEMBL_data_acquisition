"""Minimal CLI for generating dummy activity data."""

from __future__ import annotations

import argparse
import logging
from collections.abc import Sequence

from library.activities import get_activities

logger = logging.getLogger(__name__)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments.

    Parameters
    ----------
    argv
        Optional sequence of arguments for unit testing.

    Returns
    -------
    argparse.Namespace
        Parsed arguments containing ``limit`` and ``dry_run``.
    """
    parser = argparse.ArgumentParser(description="Generate dummy activity data")
    parser.add_argument(
        "--limit",
        type=int,
        default=10,
        help="Maximum number of activity rows to emit",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Log actions without generating data",
    )
    return parser.parse_args(argv)


def run(limit: int, dry_run: bool) -> int:
    """Execute the activity generation pipeline.

    Parameters
    ----------
    limit
        Number of activities to generate.
    dry_run
        If ``True``, only log the intention without generating data.

    Returns
    -------
    int
        Zero on success, non-zero on failure.
    """
    if dry_run:
        logger.info("dry run: would generate %d activities", limit)
        return 0

    activities = get_activities(limit)
    logger.info("generated %d activities", len(activities))
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry point."""
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(message)s")
    return run(limit=args.limit, dry_run=args.dry_run)


if __name__ == "__main__":
    raise SystemExit(main())

"""Minimal CLI for generating dummy activity data."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

# Allow running the script directly via ``python scripts/get_activities.py``
# by ensuring the repository root is on ``sys.path`` when the module is executed
# outside of the ``scripts`` package. This mirrors the behaviour of installing
# the project in editable mode.
if __package__ in {None, ""}:
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from library import cli
from library.activities import get_activities
from library.log import logger


def parse_args(
    argv: Sequence[str] | None = None,
) -> tuple[argparse.Namespace, cli.LoggerConfig]:
    """Parse command-line arguments.

    Parameters
    ----------
    argv
        Optional sequence of arguments for unit testing.

    Returns
    -------
    tuple[argparse.Namespace, cli.LoggerConfig]
        Parsed arguments and logging configuration.
    """
    parser, log_cfg = cli.build_parser("Generate dummy activity data", column="id")
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
    return parser.parse_args(argv), log_cfg


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
    args, log_cfg = parse_args(argv)
    log_cfg.level = args.log_level
    cli.configure_logger(log_cfg)
    return run(limit=args.limit, dry_run=args.dry_run)


if __name__ == "__main__":
    raise SystemExit(main())

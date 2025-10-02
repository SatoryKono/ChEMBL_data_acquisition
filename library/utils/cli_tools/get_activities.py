"""Minimal CLI for generating dummy activity data."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

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

    def _limit(value: str) -> int:
        """Return ``value`` validated as a non-negative integer."""

        if value == "0":
            return 0
        return cli.positive_int(value)

    parser.add_argument(
        "--limit",
        type=_limit,
        default=10,
        help="Maximum number of activity rows to emit",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Log actions without generating data",
    )
    args = parser.parse_args(argv)
    input_path = getattr(args, "input_csv", None)
    output_stem = Path(input_path).stem if input_path else None
    cli.prepare_io_paths(args, output_stem=output_stem)
    return args, log_cfg


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
        logger.info("dry_run", limit=limit)
        return 0

    activities = get_activities(limit)
    logger.info("generated", count=len(activities))
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry point.

    Notes
    -----
    Path options ``--base-path``, ``--input-dir`` and ``--output-dir`` are
    respected when resolving CLI inputs.
    """
    args, log_cfg = parse_args(argv)
    log_cfg.level = args.log_level
    cli.configure_logger(log_cfg)
    return run(limit=args.limit, dry_run=args.dry_run)


if __name__ == "__main__":
    raise SystemExit(main())

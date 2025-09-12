"""CLI for chunked CSV copying with checkpoint resume."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from library.chunk_io import process_csv_chunks
from library.cli import LoggerConfig, add_common_arguments, configure_logger
from library.config import Config, ensure_dirs
from library.io import default_output_path
from library.log import logger


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Copy CSV files in chunks with checkpointing",
    )
    add_common_arguments(parser)
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=1000,
        help="Number of rows processed per chunk",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional maximum number of rows to process",
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=Path("checkpoint.json"),
        help="Path to checkpoint file",
    )
    return parser, LoggerConfig(level="INFO")


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the chunked copy operation."""
    try:
        if args.chunk_size <= 0:
            raise SystemExit("--chunk-size must be a positive integer")
        if args.limit is not None and args.limit <= 0:
            raise SystemExit("--limit must be a positive integer")

        output = args.output_csv or default_output_path(args.input_csv, cfg.io)
        ensure_dirs(cfg)
        rows = process_csv_chunks(
            args.input_csv,
            output,
            cfg=cfg.io,
            chunk_size=args.chunk_size,
            limit=args.limit,
            checkpoint_path=args.checkpoint,
            sep=args.sep,
            encoding=args.encoding,
        )
        logger.info("rows_processed", extra={"rows": rows})
        return 0
    except Exception as exc:  # pragma: no cover - defensive
        logger.exception("run_fail", exc=exc)
        return 1


def main(argv: Sequence[str] | None = None) -> int:
    """Run the chunked copy command-line interface.

    Parameters
    ----------
    argv : Sequence[str] | None, optional
        Command-line arguments. If ``None``, values are read from :data:`sys.argv`.

    Returns
    -------
    int
        ``0`` for success, ``1`` if an error occurred during processing.
    """

    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    log_cfg.level = args.log_level
    configure_logger(log_cfg)
    cfg = Config()
    return run(cfg, args)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

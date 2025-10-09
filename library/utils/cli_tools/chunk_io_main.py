"""CLI for chunked CSV copying with checkpoint resume."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

from library import cli
from library.common.chunk_io import process_csv_chunks
from library.cli import (
    LoggerConfig,
    add_common_arguments,
    configure_logger,
    create_logger_config,
    path_argument,
)
from library.config import (
    Config,
    ConfigError,
    DEFAULT_CONFIG_PATH,
    ensure_dirs,
    print_config,
)
from library.io.paths import default_output_path
from library.common.log import logger


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Copy CSV files in chunks with checkpointing",
    )
    add_common_arguments(parser)
    parser.add_argument(
        "--chunk-size",
        type=cli.positive_int,
        default=1000,
        help="Number of rows processed per chunk",
    )
    parser.add_argument(
        "--limit",
        type=cli.positive_int,
        default=None,
        help="Optional maximum number of rows to process",
    )
    parser.add_argument(
        "--checkpoint",
        type=path_argument,
        default=Path("checkpoint.json"),
        help="Path to checkpoint file",
    )
    parser.add_argument(
        "--config",
        dest="config",
        type=path_argument,
        default=DEFAULT_CONFIG_PATH,
        help=f"YAML configuration file (default: {DEFAULT_CONFIG_PATH})",
    )
    parser.add_argument(
        "--print-config",
        action="store_true",
        help="Print effective configuration and exit",
    )
    run_id_default = parser.get_default("run_id")
    if run_id_default in (None, argparse.SUPPRESS):
        run_id_value = None
    else:
        run_id_value = str(run_id_default)
    return parser, create_logger_config("INFO", run_id=run_id_value)


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the chunked copy operation."""
    try:
        output = Path(args.output_csv or default_output_path(args.input_csv, cfg.io))
        parent = output.parent
        if not parent.exists():
            if cfg.io.exist_ok:
                parent.mkdir(parents=True, exist_ok=True)
            else:
                logger.error(
                    "output_directory_missing",
                    directory=str(parent),
                    output=str(output),
                )
                return 1
        elif not parent.is_dir():
            logger.error(
                "output_directory_not_directory",
                directory=str(parent),
                output=str(output),
            )
            return 1
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
            ensure_directory=False,
        )
        logger.info("rows_processed", rows=rows)
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

    Notes
    -----
    The command honours ``--base-path``, ``--input-dir`` and ``--output-dir``
    for resolving relative file locations.
    """

    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    input_path = getattr(args, "input_csv", None)
    output_stem = Path(input_path).stem if input_path else None
    cli.prepare_io_paths(args, output_stem=output_stem)
    run_id_value = getattr(args, "run_id", None)
    if isinstance(run_id_value, str):
        run_id_value = run_id_value.strip() or None
    if run_id_value is not None:
        log_cfg.run_id = run_id_value
    log_cfg.level = args.log_level
    configure_logger(log_cfg)
    try:
        cfg = cli.apply_config_overrides(args, parser, args.config)
    except (ConfigError, FileNotFoundError, ValueError) as exc:
        logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        return 1
    if args.print_config:
        print_config(cfg)
        return 0
    return run(cfg, args)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

"""CLI entry point for table quality analysis."""

from __future__ import annotations

import argparse
from argparse import BooleanOptionalAction
from collections.abc import Sequence
from pathlib import Path

from library import cli, io
from library.cli import (
    LoggerConfig,
    configure_logger,
)
from library.cli import (
    build_parser as base_parser,
)
from library.cli_utils import ensure_run_id
from library.common.log import logger
from library.config import Config, ConfigError, ensure_dirs, print_config
from library.table_quality import analyze_table_quality


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the profiling workflow.

    Parameters
    ----------
    cfg : Config
        Application configuration.
    args:
        Parsed command-line arguments.

    Returns
    -------
    int
        ``0`` on success, non-zero on failure.

    """
    try:
        try:
            df = io.read_csv(
                args.input_csv,
                cfg=cfg.io,
                sep=args.sep,
                encoding=args.encoding,
                dtype=str,
            )
        except io.CsvReadError as exc:
            logger.error(
                "input_read_failed",
                error=str(exc.original_error),
                path=str(args.input_csv),
                encoding=args.encoding,
                sep=args.sep,
            )
            return 1
        except FileNotFoundError as exc:
            logger.error(
                "input_read_failed",
                error=str(exc),
                path=str(args.input_csv),
                encoding=args.encoding,
                sep=args.sep,
            )
            return 1

        doc_cfg = cfg.system.doc_quality
        enable = getattr(args, "doc_quality_enable", None)
        if enable is None:
            enable = doc_cfg.enable
        if not enable:
            logger.info("doc_quality_disabled", table_name=args.table_name)
            return 0

        sample_rows = getattr(args, "sample_rows", None)
        if sample_rows is None:
            sample_rows = doc_cfg.sample_rows
        include_columns = getattr(args, "include_columns", None)
        if include_columns is None:
            include_columns = doc_cfg.include_columns
        exclude_columns = getattr(args, "exclude_columns", None)
        if exclude_columns is None:
            exclude_columns = doc_cfg.exclude_columns

        output_dir = args.output_csv

        if output_dir.suffix:
            logger.error(
                "output_directory_invalid",
                path=str(output_dir),
                reason="has_suffix",
            )
            return 1

        if output_dir.exists():
            if not output_dir.is_dir():
                logger.error(
                    "output_directory_invalid",
                    path=str(output_dir),
                    reason="not_directory",
                )
                return 1
        else:
            if not cfg.io.exist_ok:
                logger.error(
                    "output_directory_missing",
                    path=str(output_dir),
                    exist_ok=cfg.io.exist_ok,
                )
                return 1
            output_dir.mkdir(parents=True, exist_ok=True)

        destination_dir = output_dir.resolve()
        analyze_table_quality(
            df,
            table_name=args.table_name,
            destination_dir=destination_dir,
            sample_rows=sample_rows,
            include_columns=include_columns,
            exclude_columns=exclude_columns,
        )
        return 0
    except Exception as exc:  # pragma: no cover - defensive
        logger.exception("run_fail", exc=exc)
        return 1


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser."""
    parser, log_cfg = base_parser(
        "Table quality analysis (reads CSV data as strings to preserve identifiers)",
        column="chembl_id",
    )
    parser.add_argument(
        "--table-name",
        required=True,
        help="Base name used for output report files",
    )
    parser.add_argument(
        "--doc-quality-enable",
        dest="doc_quality_enable",
        action=BooleanOptionalAction,
        default=None,
        help="Enable or disable table profiling (config: system.doc_quality.enable)",
    )
    parser.add_argument(
        "--sample-rows",
        dest="sample_rows",
        type=cli.positive_int,
        default=None,
        help="Limit profiling to the first N rows",
    )
    parser.add_argument(
        "--include-columns",
        nargs="+",
        default=None,
        help="Only analyse the specified columns",
    )
    parser.add_argument(
        "--exclude-columns",
        nargs="+",
        default=None,
        help="Exclude the specified columns from analysis",
    )
    parser.set_defaults(func=run, output_csv=Path("."), encoding="utf-8-sig")

    for action in parser._actions:
        if action.dest == "output_csv":
            action.help = "Directory for profiling reports (must exist when cfg.io.exist_ok is false)"
            action.metavar = "OUTPUT_DIR"
            break
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults.

    Notes
    -----
    Relative paths honour ``--base-path``, ``--input-dir`` and ``--output-dir``.
    """
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    input_path = getattr(args, "input_csv", None)
    output_stem = Path(input_path).stem if input_path else None
    cli.prepare_io_paths(args, output_stem=output_stem)
    ensure_run_id(args, parser, log_cfg)
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline_start", run_id=log_cfg.run_id)
    try:
        cfg: Config = cli.apply_config_overrides(
            args,
            parser,
            args.config,
            mapping={
                "doc_quality_enable": "system.doc_quality.enable",
                "sample_rows": "system.doc_quality.sample_rows",
                "include_columns": "system.doc_quality.include_columns",
                "exclude_columns": "system.doc_quality.exclude_columns",
            },
        )
    except (ConfigError, FileNotFoundError, ValueError) as exc:
        logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1

    try:
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        logger = configure_logger(log_cfg)
    except (ValueError, TypeError) as exc:
        logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error(
            "directory_setup_failed",
            error=str(exc),
            output=str(args.output_csv),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    exit_code: int = args.func(cfg, args)
    if exit_code == 0:
        logger.info("pipeline_done", run_id=log_cfg.run_id)
    else:
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

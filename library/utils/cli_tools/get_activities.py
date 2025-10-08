"""Minimal CLI for generating dummy activity data."""

from __future__ import annotations

import argparse
from collections.abc import Iterable, Sequence
from pathlib import Path

import pandas as pd

from library import cli, io
from library.common.log import logger
from library.config import (
    Config,
    ConfigError,
    DEFAULT_CONFIG_PATH,
    ensure_dirs,
    print_config,
)
from library.pipelines.activity import get_activities


def parse_args(
    argv: Sequence[str] | None = None,
) -> tuple[argparse.ArgumentParser, argparse.Namespace, cli.LoggerConfig]:
    """Return parser, parsed arguments and logging configuration."""

    parser = argparse.ArgumentParser(description="Generate dummy activity data")
    cli.add_common_arguments(parser)
    parser.add_argument(
        "--config",
        dest="config",
        type=cli.path_argument,
        default=DEFAULT_CONFIG_PATH,
        help=f"YAML configuration file (default: {DEFAULT_CONFIG_PATH})",
    )
    parser.add_argument(
        "--print-config",
        action="store_true",
        help="Print effective configuration and exit",
    )
    log_cfg = cli.create_logger_config(parser.get_default("log_level"))

    def _limit(value: str) -> int:
        """Return ``value`` validated as a non-negative integer."""

        try:
            parsed = int(value)
        except ValueError as exc:  # pragma: no cover - handled by argparse
            raise argparse.ArgumentTypeError(
                "limit must be a non-negative integer"
            ) from exc
        if parsed < 0:
            raise argparse.ArgumentTypeError("limit must be a non-negative integer")
        return parsed

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
    return parser, args, log_cfg


def _frame_from_records(records: Iterable[dict[str, object]]) -> pd.DataFrame:
    """Return a dataframe built from synthetic activity records."""

    frame = pd.DataFrame(list(records))
    if frame.empty:
        return pd.DataFrame(columns=["activity_id"], dtype="Int64")
    return frame


def _write_output(frame: pd.DataFrame, output_path: Path, *, cfg: Config) -> Path:
    """Persist ``frame`` and accompanying metadata to ``output_path``."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    written = io.write_csv(frame, output_path, cfg=cfg)
    io.write_meta_yaml(
        written,
        cfg=cfg,
        columns=list(frame.columns),
        dtypes={col: str(dtype) for col, dtype in frame.dtypes.items()},
    )
    return written


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the activity generation pipeline."""

    if args.dry_run:
        logger.info("dry_run", limit=args.limit)
        return 0

    try:
        frame = _frame_from_records(get_activities(args.limit))
    except ValueError as exc:
        logger.error("invalid_arguments", error=str(exc))
        return 1

    output_candidate = getattr(args, "output_csv", None)
    input_candidate = getattr(args, "input_csv", None)
    if output_candidate in (None, argparse.SUPPRESS):
        fallback_input = input_candidate or Path("activities.csv")
        output_path = io.default_output_path(fallback_input, cfg.io)
        args.output_csv = output_path
        setattr(args, "final_out", output_path)
    else:
        output_path = Path(output_candidate)

    written_path = _write_output(frame, output_path, cfg=cfg)
    logger.info(
        "activity_generated",
        count=int(frame.shape[0]),
        output=str(written_path),
    )
    logger.info(
        "generated",
        output=str(written_path),
        count=int(frame.shape[0]),
    )
    logger.info(
        "activity_pipeline_done",
        output=str(written_path),
        rows=int(frame.shape[0]),
    )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry point."""

    parser, args, log_cfg = parse_args(argv)
    log_cfg.level = args.log_level
    structured_logger = cli.configure_logger(log_cfg)
    if structured_logger is None:  # pragma: no cover - exercised via CLI tests
        structured_logger = logger
    structured_logger.info("pipeline_start", run_id=log_cfg.run_id)

    try:
        cfg: Config = cli.apply_config_overrides(args, parser, args.config)
    except (ConfigError, FileNotFoundError, ValueError) as exc:
        logger.error("config_error", error=str(exc), config=str(args.config))
        structured_logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1

    if args.print_config:
        print_config(cfg)
        cli.configure_logger(log_cfg)
        structured_logger.info("pipeline_done", run_id=log_cfg.run_id)
        return 0

    try:
        ensure_dirs(cfg)
    except (ValueError, TypeError) as exc:
        logger.error("config_error", error=str(exc), config=str(args.config))
        structured_logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error("setup_fail", error=str(exc))
        structured_logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1

    try:
        exit_code = run(cfg, args)
    except Exception as exc:  # pragma: no cover - defensive
        logger.exception("run_fail", exc=exc)
        structured_logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1

    if exit_code == 0:
        structured_logger.info("pipeline_done", run_id=log_cfg.run_id)
    else:
        structured_logger.info("pipeline_fail", run_id=log_cfg.run_id)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

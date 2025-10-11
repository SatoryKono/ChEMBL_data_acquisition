"""Minimal CLI for generating dummy activity data."""

from __future__ import annotations

import argparse
from collections.abc import Iterable, Sequence
from pathlib import Path
from uuid import uuid4

import pandas as pd
import yaml

from library import cli, io
from library.common.log import logger
from library.config import (
    DEFAULT_CONFIG_PATH,
    Config,
    ConfigError,
    ensure_dirs,
    print_config,
)
from library.pipelines.activity import get_activities
from library.utils.atomic import robust_replace

DEFAULT_LIMIT = 25


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
    run_id_default = parser.get_default("run_id")
    if run_id_default in (None, argparse.SUPPRESS):
        run_id_value = None
    else:
        run_id_value = str(run_id_default)
    log_cfg = cli.create_logger_config(
        parser.get_default("log_level"),
        run_id=run_id_value,
    )

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
        default=None,
        help=f"Maximum number of activity rows to emit (default: {DEFAULT_LIMIT})",
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
    run_id_value = getattr(args, "run_id", None)
    if isinstance(run_id_value, str):
        run_id_value = run_id_value.strip() or None
    if run_id_value is not None:
        log_cfg.run_id = run_id_value
    return parser, args, log_cfg


def _frame_from_records(records: Iterable[dict[str, object]]) -> pd.DataFrame:
    """Return a dataframe built from synthetic activity records."""

    frame = pd.DataFrame(list(records))
    if frame.empty:
        return pd.DataFrame(columns=["activity_id"], dtype="Int64")
    return frame


def _write_output(frame: pd.DataFrame, output_path: Path, *, cfg: Config) -> Path:
    """Persist ``frame`` and accompanying metadata to ``output_path`` atomically."""

    output_path.parent.mkdir(parents=True, exist_ok=True)

    tmp_path = output_path.with_name(f".{output_path.name}.{uuid4().hex}.tmp")
    written = io.write_csv(frame, tmp_path, cfg=cfg)
    tmp_meta = Path(f"{written}.meta.yaml")

    if not tmp_meta.exists():
        written.unlink(missing_ok=True)
        msg = f"metadata sidecar missing for temporary output: {tmp_meta}"
        raise FileNotFoundError(msg)

    metadata = yaml.safe_load(tmp_meta.read_text(encoding="utf-8")) or {}
    expected_columns = list(frame.columns)
    expected_dtypes = {column: str(dtype) for column, dtype in frame.dtypes.items()}

    meta_columns = metadata.get("columns")
    meta_dtypes = metadata.get("dtypes")
    if meta_columns != expected_columns or meta_dtypes != expected_dtypes:
        written.unlink(missing_ok=True)
        tmp_meta.unlink(missing_ok=True)
        msg = "metadata sidecar schema mismatch for generated CSV"
        raise ValueError(msg)

    target_meta = Path(f"{output_path}.meta.yaml")

    try:
        robust_replace(written, output_path)
        robust_replace(tmp_meta, target_meta)
    except Exception:
        written.unlink(missing_ok=True)
        tmp_meta.unlink(missing_ok=True)
        raise

    # Clean up any lingering temporary metadata artefacts from the staging write.
    for orphan in output_path.parent.glob(
        f".{output_path.name}.*.tmp.meta.yaml*"
    ):
        orphan.unlink(missing_ok=True)

    return output_path


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the activity generation pipeline."""

    limit = args.limit
    if limit is None:
        limit = cfg.activity.limit if cfg.activity.limit is not None else DEFAULT_LIMIT
        if limit < 0:
            logger.error(
                "config_error",
                error="activity.limit must be zero or a positive integer",
                limit=limit,
            )
            return 1

    args.limit = limit

    if args.dry_run:
        logger.info("dry_run", limit=limit)
        return 0

    output_candidate = getattr(args, "output_csv", None)
    input_candidate = getattr(args, "input_csv", None)
    if output_candidate in (None, argparse.SUPPRESS):
        fallback_input = input_candidate or Path("activities.csv")
        output_path = io.default_output_path(
            fallback_input,
            cfg.io,
            date=getattr(args, "date", None),
        )
        args.output_csv = output_path
        args.final_out = output_path
    else:
        output_path = Path(output_candidate)

    if (
        output_path.exists()
        and getattr(args, "skip_existing", False)
        and not getattr(args, "force", False)
    ):
        logger.info(
            "pipeline_skip_existing", output_postprocessed=str(output_path)
        )
        return 0

    try:
        frame = _frame_from_records(get_activities(limit))
    except ValueError as exc:
        logger.error("invalid_arguments", error=str(exc))
        return 1

    written_path = _write_output(frame, output_path, cfg=cfg)
    logger.info(
        "activity_generated",
        count=int(frame.shape[0]),
        output_postprocessed=str(written_path),
    )
    logger.info(
        "generated",
        output_postprocessed=str(written_path),
        count=int(frame.shape[0]),
    )
    logger.info(
        "activity_pipeline_done",
        output_postprocessed=str(written_path),
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

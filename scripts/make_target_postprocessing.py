"""Command line interface for postprocessing target exports."""

from __future__ import annotations

if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli

    bootstrap_cli(__package__, __file__)
    package_name = "scripts"
else:  # pragma: no cover - executed when imported as a package module
    from ._bootstrap import bootstrap_cli

    bootstrap_cli(__package__, __file__)
    package_name = __package__

del bootstrap_cli

import argparse
import os
from importlib import import_module
from pathlib import Path
from typing import Sequence

import pandas as pd

from library import io  # noqa: F401 - imported for CLI parity with existing scripts
from library.cli import configure_logger, create_logger_config
from library.cli.logging import setup_cli_logging
from library.cli.parser import path_argument
from library.common.log import logger
from library.postprocess.common.types import SchemaValidationError, StepError
from library.postprocess.targets.schema import TARGET_SCHEMA, validate_targets
from library.postprocess.targets.steps import run_target_pipeline

_postprocess_common = import_module(f"{package_name}._postprocess_common")

CsvRuntimeConfig = _postprocess_common.CsvRuntimeConfig
DEFAULT_LOG_DIR = _postprocess_common.DEFAULT_LOG_DIR
LOG_DIR_ENV = _postprocess_common.LOG_DIR_ENV
event_prefix = _postprocess_common.event_prefix
export_postprocess_frame = _postprocess_common.export_postprocess_frame
generate_metrics_report = _postprocess_common.generate_metrics_report
get_csv_runtime_config = _postprocess_common.get_csv_runtime_config
get_default_log_level = _postprocess_common.get_default_log_level
get_pipeline_config = _postprocess_common.get_pipeline_config
load_input_frame = _postprocess_common.load_input_frame
run_postprocess_steps = _postprocess_common.run_postprocess_steps
validate_postprocess_frame = _postprocess_common.validate_postprocess_frame
del _postprocess_common, package_name


TABLE_NAME = "targets"
PROGRAM_NAME = Path(__file__).with_suffix("").name


def build_parser() -> argparse.ArgumentParser:
    """Create the command line parser for the target postprocessor."""

    parser = argparse.ArgumentParser(
        description="Run the target postprocessing pipeline on a CSV export.",
    )
    parser.add_argument(
        "--input",
        dest="input",
        type=path_argument,
        required=True,
        help="Input CSV file containing raw target records.",
    )
    parser.add_argument(
        "--output",
        dest="output",
        type=path_argument,
        required=True,
        help="Destination CSV file for the postprocessed targets.",
    )
    parser.add_argument(
        "--config",
        dest="config",
        type=path_argument,
        default=None,
        help="Optional YAML configuration overriding the default pipeline definition.",
    )
    parser.add_argument(
        "--log-level",
        dest="log_level",
        default=None,
        help="Logging verbosity (defaults to the pipeline configuration).",
    )
    return parser


def load_output_data(path: Path, csv_cfg: CsvRuntimeConfig) -> pd.DataFrame:
    """Load the raw target output frame."""

    return load_input_frame(TABLE_NAME, path, csv_cfg, logger=logger)

def run(args: argparse.Namespace) -> int:
    """Execute the target postprocessing pipeline using ``args``."""

    pipeline_config = getattr(args, "_pipeline_config", None)
    if pipeline_config is None:
        pipeline_config = get_pipeline_config(TABLE_NAME, getattr(args, "config", None))
    csv_cfg = getattr(args, "_csv_runtime_config", None)
    if csv_cfg is None:
        csv_cfg = get_csv_runtime_config(pipeline_config)

    input_path = Path(args.input)
    output_path = Path(args.output)
    prefix = event_prefix(TABLE_NAME)

    if not input_path.exists():
        logger.error(f"{prefix}_input_missing", input=str(input_path))
        return 1

    metrics = None
    processed: pd.DataFrame | None = None

    try:
        frame = load_output_data(input_path, csv_cfg)
        processed, metrics = run_postprocess_steps(
            TABLE_NAME,
            frame,
            run_target_pipeline,
            pipeline_config.pipeline_version,
            logger=logger,
        )
        effective_version = (
            metrics.pipeline_version if metrics and metrics.pipeline_version else pipeline_config.pipeline_version
        )
        validated = validate_postprocess_frame(
            TABLE_NAME,
            processed,
            validate_targets,
            TARGET_SCHEMA,
            effective_version,
            logger=logger,
        )
        export_postprocess_frame(
            TABLE_NAME,
            validated,
            output_path,
            csv_cfg,
            TARGET_SCHEMA,
            logger=logger,
        )
        processed = validated
    except (SchemaValidationError, StepError) as exc:
        logger.exception(f"{prefix}_failed", exc=exc)
        return 1
    except Exception as exc:  # pragma: no cover - defensive guard
        logger.exception(f"{prefix}_unexpected_error", exc=exc)
        return 1

    if metrics is not None:
        summary = metrics.summary()
        summary["output"] = str(output_path)
        logger.info(f"{prefix}_summary", **summary)

    extras = {"input": str(input_path), "output": str(output_path)}
    effective_version = (
        metrics.pipeline_version if metrics and metrics.pipeline_version else pipeline_config.pipeline_version
    )
    generate_metrics_report(
        TABLE_NAME,
        output_path,
        csv_cfg,
        run_target_pipeline,
        pipeline_version=effective_version,
        extras=extras,
        logger=logger,
        pipeline_metrics=metrics,
    )

    output_rows = (
        int(metrics.output_rows)
        if metrics and metrics.output_rows is not None
        else int(processed.shape[0])
        if processed is not None
        else 0
    )
    output_columns = (
        int(metrics.output_columns)
        if metrics and metrics.output_columns is not None
        else int(processed.shape[1])
        if processed is not None
        else 0
    )
    logger.info(
        f"{prefix}_done",
        output=str(output_path),
        rows=output_rows,
        columns=output_columns,
    )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point wiring argument parsing and logging configuration."""

    parser = build_parser()
    args = parser.parse_args(argv)

    pipeline_config = get_pipeline_config(TABLE_NAME, getattr(args, "config", None))
    csv_cfg = get_csv_runtime_config(pipeline_config)
    log_level = (args.log_level or get_default_log_level(pipeline_config)).upper()
    log_cfg = create_logger_config(log_level)

    log_dir_value = os.environ.get(LOG_DIR_ENV)
    if log_dir_value:
        log_dir = Path(log_dir_value)
    else:
        log_dir = DEFAULT_LOG_DIR

    with setup_cli_logging(PROGRAM_NAME, log_cfg, date_str=None, log_dir=log_dir) as logging_ctx:
        configure_logger(logging_ctx.log_cfg)
        setattr(args, "_pipeline_config", pipeline_config)
        setattr(args, "_csv_runtime_config", csv_cfg)
        setattr(args, "_log_level", log_level)
        exit_code = run(args)

    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())


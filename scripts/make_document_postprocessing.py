"""Command line interface for postprocessing document exports."""

from __future__ import annotations

if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli
else:  # pragma: no cover - executed when imported as a package module
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli

import argparse
import os
from pathlib import Path
from typing import Sequence

import pandas as pd

from library import io  # noqa: F401 - imported for CLI parity with existing scripts
from library.cli import configure_logger, create_logger_config
from library.cli.logging import setup_cli_logging
from library.cli.parser import path_argument
from library.common.log import logger
from library.postprocess.common.logging import PipelineRunMetrics
from library.postprocess.common.types import SchemaValidationError, StepError
from library.postprocess.documents.schema import DOCUMENT_SCHEMA, validate_documents
from library.postprocess.documents.steps import run_document_pipeline

try:
    from _postprocess_common import (
        CsvRuntimeConfig,
        DEFAULT_LOG_DIR,
        LOG_DIR_ENV,
        export_postprocess_frame,
        generate_metrics_report,
        get_csv_runtime_config,
        get_default_log_level,
        get_pipeline_config,
        load_input_frame,
        run_postprocess_steps,
        validate_postprocess_frame,
    )
except ModuleNotFoundError:  # pragma: no cover - fallback when executed as a package module
    from ._postprocess_common import (
        CsvRuntimeConfig,
        DEFAULT_LOG_DIR,
        LOG_DIR_ENV,
        export_postprocess_frame,
        generate_metrics_report,
        get_csv_runtime_config,
        get_default_log_level,
        get_pipeline_config,
        load_input_frame,
        run_postprocess_steps,
        validate_postprocess_frame,
    )


TABLE_NAME = "documents"
PROGRAM_NAME = Path(__file__).with_suffix("").name


def build_parser() -> argparse.ArgumentParser:
    """Create the command line parser for the document postprocessor."""

    parser = argparse.ArgumentParser(
        description="Run the document postprocessing pipeline on a CSV export.",
    )
    parser.add_argument(
        "--input",
        dest="input",
        type=path_argument,
        required=True,
        help="Input CSV file containing raw document records.",
    )
    parser.add_argument(
        "--output",
        dest="output",
        type=path_argument,
        required=True,
        help="Destination CSV file for the postprocessed documents.",
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
    """Load the raw document output frame."""

    return load_input_frame(TABLE_NAME, path, csv_cfg, logger=logger)


def apply_postprocessing_steps(
    df: pd.DataFrame,
    *,
    pipeline_version: str | None,
) -> tuple[pd.DataFrame, PipelineRunMetrics]:
    """Execute the configured document postprocessing pipeline."""

    return run_postprocess_steps(
        TABLE_NAME,
        df,
        run_document_pipeline,
        pipeline_version,
        logger=logger,
    )


def validate_output_schema(
    df: pd.DataFrame,
    *,
    pipeline_version: str | None,
) -> pd.DataFrame:
    """Validate and reorder the DataFrame according to the document schema."""

    return validate_postprocess_frame(
        TABLE_NAME,
        df,
        validate_documents,
        DOCUMENT_SCHEMA,
        pipeline_version,
        logger=logger,
    )


def save_output_data(
    df: pd.DataFrame,
    output_path: Path,
    csv_cfg: CsvRuntimeConfig,
) -> Path:
    """Persist ``df`` deterministically to ``output_path``."""

    return export_postprocess_frame(
        TABLE_NAME,
        df,
        output_path,
        csv_cfg,
        DOCUMENT_SCHEMA,
        logger=logger,
    )


def run(args: argparse.Namespace) -> int:
    """Execute the document postprocessing pipeline using ``args``."""

    pipeline_config = getattr(args, "_pipeline_config", None)
    if pipeline_config is None:
        pipeline_config = get_pipeline_config(TABLE_NAME, getattr(args, "config", None))
    csv_cfg = getattr(args, "_csv_runtime_config", None)
    if csv_cfg is None:
        csv_cfg = get_csv_runtime_config(pipeline_config)

    input_path = Path(args.input)
    output_path = Path(args.output)
    event_prefix = f"{TABLE_NAME}_postprocess"

    if not input_path.exists():
        logger.error(f"{event_prefix}_input_missing", input=str(input_path))
        return 1

    metrics: PipelineRunMetrics | None = None

    try:
        frame = load_output_data(input_path, csv_cfg)
        processed, metrics = apply_postprocessing_steps(
            frame,
            pipeline_version=pipeline_config.pipeline_version,
        )
        validated = validate_output_schema(
            processed,
            pipeline_version=metrics.pipeline_version,
        )
        save_output_data(validated, output_path, csv_cfg)
    except (SchemaValidationError, StepError) as exc:
        logger.exception(f"{event_prefix}_failed", exc=exc)
        return 1
    except Exception as exc:  # pragma: no cover - defensive guard
        logger.exception(f"{event_prefix}_unexpected_error", exc=exc)
        return 1

    if metrics is not None:
        summary = metrics.summary()
        summary["output"] = str(output_path)
        logger.info(f"{event_prefix}_summary", **summary)

    extras = {"input": str(input_path), "output": str(output_path)}
    generate_metrics_report(
        TABLE_NAME,
        output_path,
        csv_cfg,
        run_document_pipeline,
        pipeline_version=metrics.pipeline_version if metrics else None,
        extras=extras,
        logger=logger,
    )

    logger.info(
        f"{event_prefix}_done",
        output=str(output_path),
        rows=int(validated.shape[0]),
        columns=int(validated.shape[1]),
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


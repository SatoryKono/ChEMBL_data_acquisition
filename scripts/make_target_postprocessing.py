"""Command line interface for postprocessing target exports."""

from __future__ import annotations

# ruff: noqa: E402  # bootstrap alters import order for script compatibility
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
from collections.abc import Mapping, Sequence
from importlib import import_module
from pathlib import Path
from uuid import NAMESPACE_URL, uuid5

import pandas as pd

from library import io  # noqa: F401 - imported for CLI parity with existing scripts
from library.cli import configure_logger, create_logger_config
from library.cli.logging import setup_cli_logging
from library.cli.parser import path_argument
from library.cli_utils import resolve_invocation
from library.common.log import logger
from library.pipelines.common.metadata import get_pipeline_version
from library.postprocessing.common.config import (
    PipelineConfig,
    normalize_pipeline_version,
)
from library.postprocessing.common.logging import PipelineRunMetrics
from library.postprocessing.common.types import SchemaValidationError, StepError
from library.postprocessing.targets import (
    run_target_pipeline as run_target_postprocess,
)
from library.postprocessing.targets import steps as target_steps
from library.postprocessing.targets.schema import TARGET_SCHEMA, validate_targets

_postprocess_common = import_module(f"{package_name}._postprocess_common")

CsvRuntimeConfig = _postprocess_common.CsvRuntimeConfig
DEFAULT_LOG_DIR = _postprocess_common.DEFAULT_LOG_DIR
LOG_DIR_ENV = _postprocess_common.LOG_DIR_ENV
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
    parser.add_argument(
        "--run-id",
        dest="run_id",
        default=os.environ.get("CHEMBL_DA_RUN_ID"),
        help="Override the run identifier used for logging",
    )
    return parser


def load_output_data(path: Path, csv_cfg: CsvRuntimeConfig) -> pd.DataFrame:
    """Load the raw target output frame."""

    return load_input_frame(TABLE_NAME, path, csv_cfg, logger=logger)


def apply_postprocessing_steps(
    df: pd.DataFrame,
    *,
    pipeline_version: str | None,
) -> tuple[pd.DataFrame, PipelineRunMetrics]:
    """Execute the configured target postprocessing pipeline."""

    return run_postprocess_steps(
        TABLE_NAME,
        df,
        run_target_postprocess,
        pipeline_version,
        logger=logger,
    )


def validate_output_schema(
    df: pd.DataFrame,
    *,
    pipeline_version: str | None,
) -> pd.DataFrame:
    """Validate and reorder the DataFrame according to the target schema."""

    return validate_postprocess_frame(
        TABLE_NAME,
        df,
        validate_targets,
        TARGET_SCHEMA,
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
        TARGET_SCHEMA,
        logger=logger,
    )


def resolve_pipeline_version(
    pipeline_config: PipelineConfig,
    *,
    override: str | None = None,
) -> str:
    """Return the effective pipeline version for the current execution."""

    candidate = normalize_pipeline_version(override)
    if candidate is not None:
        return candidate

    candidate = normalize_pipeline_version(pipeline_config.pipeline_version)
    if candidate is not None:
        return candidate

    fallback = _pipeline_version_from_defaults(pipeline_config.params)
    candidate = normalize_pipeline_version(fallback)
    if candidate is not None:
        return candidate

    return get_pipeline_version()


def _pipeline_version_from_defaults(params: Mapping[str, object] | None) -> str | None:
    """Return the pipeline version declared under ``params.defaults`` when present."""

    if not params:
        return None

    params_map = dict(params)
    defaults = params_map.get("defaults")
    if not isinstance(defaults, Mapping):
        return None

    defaults_map = dict(defaults)
    value = defaults_map.get("pipeline_version")
    if value is None:
        return None
    return str(value)


def run(args: argparse.Namespace) -> int:
    """Execute the target postprocessing pipeline using ``args``."""

    pipeline_config = getattr(args, "_pipeline_config", None)
    if pipeline_config is None:
        pipeline_config = get_pipeline_config(TABLE_NAME, getattr(args, "config", None))
    csv_cfg = getattr(args, "_csv_runtime_config", None)
    if csv_cfg is None:
        csv_cfg = get_csv_runtime_config(pipeline_config)

    target_steps.PIPELINE_CONFIG = pipeline_config
    target_steps.PIPELINE_STEPS = pipeline_config.step_definitions()
    resolved_pipeline_version = resolve_pipeline_version(pipeline_config)

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
            pipeline_version=resolved_pipeline_version,
        )
        effective_version = metrics.pipeline_version if metrics else None
        validated = validate_output_schema(
            processed,
            pipeline_version=effective_version,
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
        run_target_postprocess,
        pipeline_version=resolved_pipeline_version,
        extras=extras,
        logger=logger,
        pipeline_metrics=metrics,
    )

    logger.info(
        f"{event_prefix}_done",
        output=str(output_path),
        rows=(
            int(metrics.output_rows)
            if metrics and metrics.output_rows is not None
            else None
        ),
        columns=(
            int(metrics.output_columns)
            if metrics and metrics.output_columns is not None
            else None
        ),
    )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point wiring argument parsing and logging configuration."""

    parser = build_parser()
    args = parser.parse_args(argv)

    pipeline_config = get_pipeline_config(TABLE_NAME, getattr(args, "config", None))
    csv_cfg = get_csv_runtime_config(pipeline_config)
    log_level = (args.log_level or get_default_log_level(pipeline_config)).upper()
    invocation = resolve_invocation(parser.prog, argv)
    run_id_value = getattr(args, "run_id", None)
    if isinstance(run_id_value, str):
        run_id_value = run_id_value.strip() or None
    if not run_id_value:
        descriptor = "\n".join(
            [
                *invocation,
                f"input={Path(args.input).resolve()}",
                f"output={Path(args.output).resolve()}",
            ]
        )
        run_id_value = uuid5(NAMESPACE_URL, descriptor).hex
    log_cfg = create_logger_config(log_level, run_id=run_id_value)

    log_dir_value = os.environ.get(LOG_DIR_ENV)
    if log_dir_value:
        log_dir = Path(log_dir_value)
    else:
        log_dir = DEFAULT_LOG_DIR

    with setup_cli_logging(
        PROGRAM_NAME, log_cfg, date_str=None, log_dir=log_dir
    ) as logging_ctx:
        configure_logger(logging_ctx.log_cfg)
        args._pipeline_config = pipeline_config
        args._csv_runtime_config = csv_cfg
        args._log_level = log_level
        exit_code = run(args)

    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

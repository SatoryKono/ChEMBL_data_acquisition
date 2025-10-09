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
from datetime import datetime, timezone
from importlib import import_module
from pathlib import Path
from time import perf_counter

import pandas as pd
from pandera.errors import SchemaErrors

from library import io  # noqa: F401 - imported for CLI parity with existing scripts
from library.cli import configure_logger, create_logger_config
from library.cli.logging import setup_cli_logging
from library.cli.parser import path_argument
from library.common.csv_utils import write_csv_deterministic
from library.common.log import logger
from library.pipelines.common import add_pipeline_metadata
from library.pipelines.target import postprocessing as target_postprocessing
from library.postprocess.common.logging import (
    PipelineRunMetrics,
    ValidationMetrics,
    execute_step,
)
from library.postprocess.common.types import SchemaValidationError, StepError
from library.postprocess.target.export import prepare_targets_for_schema
from library.schemas import TargetsSchema, normalize_targets
from library.validation import validate_targets
from library.schemas.targets import TARGETS_COLUMN_ORDER

_postprocess_common = import_module(f"{package_name}._postprocess_common")

CsvRuntimeConfig = _postprocess_common.CsvRuntimeConfig
DEFAULT_LOG_DIR = _postprocess_common.DEFAULT_LOG_DIR
LOG_DIR_ENV = _postprocess_common.LOG_DIR_ENV
generate_metrics_report = _postprocess_common.generate_metrics_report
get_csv_runtime_config = _postprocess_common.get_csv_runtime_config
get_default_log_level = _postprocess_common.get_default_log_level
get_pipeline_config = _postprocess_common.get_pipeline_config
load_input_frame = _postprocess_common.load_input_frame
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


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _run_postprocess(
    df: pd.DataFrame,
    *,
    pipeline_version: str | None,
) -> tuple[pd.DataFrame, PipelineRunMetrics]:
    """Run the target postprocessing sequence and capture metrics."""

    current = df.copy(deep=True)
    pipeline_started_at = _now_iso()
    pipeline_clock = perf_counter()
    metrics = PipelineRunMetrics(
        pipeline_version=pipeline_version,
        started_at=pipeline_started_at,
        input_rows=current.shape[0],
        input_columns=current.shape[1],
    )
    missing_optional: set[str] = set()

    def _run_step(name: str, func) -> None:
        nonlocal current
        current, step_metrics = execute_step(name, func, current, logger=logger)
        metrics.steps.append(step_metrics)

    _run_step("postprocess_targets", target_postprocessing.postprocess_targets)
    _run_step("finalise_targets", target_postprocessing.finalise_targets)
    _run_step("normalize_targets", normalize_targets)

    def _add_metadata(frame: pd.DataFrame) -> pd.DataFrame:
        enriched = add_pipeline_metadata(frame)
        if pipeline_version:
            enriched["pipeline_version"] = str(pipeline_version)
        return enriched

    _run_step("add_pipeline_metadata", _add_metadata)

    def _prepare(frame: pd.DataFrame) -> pd.DataFrame:
        prepared, missing_required, missing_optional_columns = prepare_targets_for_schema(
            frame
        )
        if missing_required:
            message = ", ".join(sorted(missing_required))
            raise SchemaValidationError(
                "prepare_targets_for_schema",
                f"missing required columns: {message}",
            )
        if missing_optional_columns:
            missing_optional.update(missing_optional_columns)
        return prepared

    _run_step("prepare_targets_for_schema", _prepare)

    validation_started_at = _now_iso()
    validation_clock = perf_counter()
    try:
        validation = validate_targets(current, return_result=True)
    except SchemaErrors as exc:  # pragma: no cover - pandera raises SchemaErrors
        raise SchemaValidationError("TargetsSchema", str(exc)) from exc
    if not validation.failure_cases.empty:
        failure_count = len(validation.failure_cases)
        sample = validation.failure_cases.head(3).to_dict("records")
        message = (
            f"{failure_count} rows failed TargetsSchema validation: "
            f"{sample}"
        )
        raise SchemaValidationError("TargetsSchema", message)
    current = validation.data
    validation_duration = perf_counter() - validation_clock
    metrics.validation = ValidationMetrics(
        schema="TargetsSchema",
        started_at=validation_started_at,
        completed_at=_now_iso(),
        duration_s=validation_duration,
    )

    def _fill_placeholders(frame: pd.DataFrame) -> pd.DataFrame:
        filled = frame.fillna("-")
        return filled.replace("", "-")

    _run_step("fill_missing_values", _fill_placeholders)

    def _drop_duplicates(frame: pd.DataFrame) -> pd.DataFrame:
        if "target_chembl_id" not in frame.columns:
            raise SchemaValidationError(
                "drop_duplicates",
                "missing required column 'target_chembl_id'",
            )
        return frame.drop_duplicates(subset=["target_chembl_id"], keep="first")

    _run_step("drop_duplicates", _drop_duplicates)

    def _sort(frame: pd.DataFrame) -> pd.DataFrame:
        if "target_chembl_id" not in frame.columns:
            raise SchemaValidationError(
                "sort_targets",
                "missing required column 'target_chembl_id'",
            )
        ordered = frame.sort_values(by=["target_chembl_id"], kind="mergesort")
        return ordered.reset_index(drop=True)

    _run_step("sort_targets", _sort)

    metrics.finalize(
        output_rows=current.shape[0],
        output_columns=current.shape[1],
        duration_s=perf_counter() - pipeline_clock,
    )

    if missing_optional:
        logger.info(
            f"{TABLE_NAME}_postprocess_missing_optional",
            columns=sorted(missing_optional),
        )

    return current, metrics


def save_output_data(
    df: pd.DataFrame,
    output_path: Path,
    csv_cfg: CsvRuntimeConfig,
) -> Path:
    """Persist ``df`` deterministically to ``output_path``."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    col_order = [column for column in TARGETS_COLUMN_ORDER if column in df.columns]
    return write_csv_deterministic(
        df.copy(),
        output,
        col_order=col_order,
        key_cols=["target_chembl_id"],
        chunksize=csv_cfg.chunksize,
        sep=csv_cfg.sep,
        encoding=csv_cfg.encoding,
        cfg=None,
    )


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
    event_prefix = f"{TABLE_NAME}_postprocess"

    if not input_path.exists():
        logger.error(f"{event_prefix}_input_missing", input=str(input_path))
        return 1

    metrics: PipelineRunMetrics | None = None

    try:
        frame = load_output_data(input_path, csv_cfg)
        processed, metrics = _run_postprocess(
            frame,
            pipeline_version=pipeline_config.pipeline_version,
        )
        save_output_data(processed, output_path, csv_cfg)
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
        _run_postprocess,
        pipeline_version=metrics.pipeline_version if metrics else None,
        extras=extras,
        logger=logger,
        pipeline_metrics=metrics,
    )

    output_rows = (
        int(metrics.output_rows)
        if metrics and metrics.output_rows is not None
        else int(processed.shape[0])
    )
    output_columns = (
        int(metrics.output_columns)
        if metrics and metrics.output_columns is not None
        else int(processed.shape[1])
    )
    logger.info(
        f"{event_prefix}_done",
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


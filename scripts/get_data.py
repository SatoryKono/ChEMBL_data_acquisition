"""Orchestrate all ChEMBL data acquisition pipelines via a unified CLI interface.

This module provides a single entry point that wires the existing acquisition
scripts together. ``get_data.py`` accepts common options for locating input and
output artefacts, resolves those paths once, and then invokes the following
modules in order while forwarding the computed arguments:

1. :mod:`scripts.get_document_data`
2. :mod:`scripts.get_target_data`
3. :mod:`scripts.get_assay_data`
4. :mod:`scripts.get_testitem_data`
5. :mod:`scripts.get_activity_data`

Each pipeline receives consistent ``--config``, ``--input`` and ``--output``
values alongside shared logging options. The command line interface exposes
parameters for the base data directory, distinct input/output sub-directories,
a date prefix for generated filenames and optional execution flags. The helper
functions defined below encapsulate the argument preparation for every step so
that the pipelines can be executed programmatically from other callers as well.
"""

from __future__ import annotations

import argparse
import sys
import uuid
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Callable, Sequence

if __package__ in {None, ""}:
    _PROJECT_ROOT = Path(__file__).resolve().parent.parent
    if str(_PROJECT_ROOT) not in sys.path:
        sys.path.insert(0, str(_PROJECT_ROOT))

from scripts import (
    get_activity_data,
    get_assay_data,
    get_document_data,
    get_target_data,
    get_testitem_data,
)

from library.logging_setup import Logger, LoggerConfig, configure_logger


_LOGGER: Logger = Logger(LoggerConfig())

_DEFAULT_INPUT_FILES = {
    "document": "document.csv",
    "target": "target.csv",
    "assay": "assay.csv",
    "testitem": "testitem.csv",
    "activity": "activity.csv",
}

_DEFAULT_OUTPUT_STEMS = {
    "document": "documents",
    "target": "targets",
    "assay": "assays",
    "testitem": "testitems",
    "activity": "activities",
}


@dataclass(frozen=True)
class PipelineRunConfig:
    """Resolved configuration shared across pipeline steps."""

    base_path: Path
    input_dir: Path
    output_dir: Path
    config_path: Path
    date_prefix: str
    log_level: str
    limit: int | None
    force: bool
    skip_existing: bool

    def input_path(self, name: str) -> Path:
        """Return the fully resolved path for ``name`` in the input directory."""

        filename = _DEFAULT_INPUT_FILES[name]
        return self.input_dir / filename

    def output_path(self, name: str) -> Path:
        """Return the fully resolved path for ``name`` in the output directory."""

        stem = _DEFAULT_OUTPUT_STEMS[name]
        filename = f"output.{stem}_{self.date_prefix}.csv"
        return self.output_dir / filename


@dataclass(frozen=True)
class PipelineStep:
    """Describe a single pipeline invocation."""

    name: str
    main: Callable[[Sequence[str] | None], int]
    subcommand: str | None

    def build_arguments(
        self, cfg: PipelineRunConfig, output_path: Path | None = None
    ) -> list[str]:
        """Return CLI arguments forwarded to the wrapped ``main`` function."""

        input_csv = cfg.input_path(self.name)
        output_csv = output_path if output_path is not None else cfg.output_path(self.name)
        args = ["--config", str(cfg.config_path), "--input", str(input_csv)]
        args.extend(["--output", str(output_csv)])
        args.extend(["--log-level", cfg.log_level])
        if cfg.limit is not None:
            args.extend(["--limit", str(cfg.limit)])
        if self.subcommand is not None:
            return [self.subcommand, *args]
        return args

    def expected_output(self, cfg: PipelineRunConfig) -> Path:
        """Return the path where the pipeline will create its CSV artefact."""

        return cfg.output_path(self.name)

    def required_input(self, cfg: PipelineRunConfig) -> Path:
        """Return the CSV that the pipeline expects as input."""

        return cfg.input_path(self.name)


_PIPELINE_STEPS: tuple[PipelineStep, ...] = (
    PipelineStep("document", get_document_data.main, "all"),
    PipelineStep("target", get_target_data.main, "all"),
    PipelineStep("assay", get_assay_data.main, None),
    PipelineStep("testitem", get_testitem_data.main, None),
    PipelineStep("activity", get_activity_data.main, None),
)


def _resolve_path(base: Path, candidate: Path) -> Path:
    """Return ``candidate`` if absolute, otherwise relative to ``base``."""

    expanded = candidate.expanduser()
    if expanded.is_absolute():
        return expanded.resolve()
    return (base / expanded).resolve()


def _configure_logging(level_name: str, *, run_id: str | None = None) -> Logger:
    """Configure structured logging for the orchestration workflow."""

    normalised = level_name.upper()
    valid_levels = {"DEBUG", "INFO", "WARN", "WARNING", "ERROR"}
    if normalised not in valid_levels:
        raise ValueError(f"invalid log level: {level_name}")

    resolved_run_id = run_id or uuid.uuid4().hex

    return configure_logger(LoggerConfig(level=normalised, run_id=resolved_run_id))


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    """Create the CLI parser and return parsed arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Run the document, target, assay, test item and activity pipelines "
            "sequentially with shared configuration."
        )
    )
    parser.add_argument(
        "--base-path",
        type=Path,
        default=Path("data"),
        help="Root directory containing input and output folders",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("input"),
        help="Directory with source CSV files relative to --base-path",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output"),
        help="Destination directory relative to --base-path",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/config.yaml"),
        help="YAML configuration shared by all pipelines",
    )
    parser.add_argument(
        "--date",
        dest="date_prefix",
        default=datetime.now(UTC).strftime("%Y%m%d"),
        help="Date prefix used to build output filenames",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging verbosity for the orchestrator and child pipelines",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of identifiers processed by each pipeline",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run pipelines even if the target file already exists",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip pipeline execution when the output file is present",
    )
    return parser.parse_args(argv)


def _prepare_config(args: argparse.Namespace) -> PipelineRunConfig:
    """Validate CLI inputs and construct :class:`PipelineRunConfig`."""

    base_path = args.base_path.expanduser().resolve()
    input_dir = _resolve_path(base_path, args.input_dir)
    output_dir = _resolve_path(base_path, args.output_dir)
    config_path = args.config.expanduser().resolve()

    if args.limit is not None and args.limit <= 0:
        raise ValueError("--limit must be a positive integer")

    if not input_dir.exists():
        raise FileNotFoundError(f"input directory not found: {input_dir}")
    if not config_path.exists():
        raise FileNotFoundError(f"configuration file not found: {config_path}")
    output_dir.mkdir(parents=True, exist_ok=True)

    return PipelineRunConfig(
        base_path=base_path,
        input_dir=input_dir,
        output_dir=output_dir,
        config_path=config_path,
        date_prefix=args.date_prefix,
        log_level=args.log_level,
        limit=args.limit,
        force=args.force,
        skip_existing=args.skip_existing,
    )


@dataclass(frozen=True)
class StepExecutionResult:
    """Summarize the outcome of invoking a pipeline step."""

    exit_code: int
    executed: bool


def _temporary_output_path(output_path: Path) -> Path:
    """Return the path used for intermediate artefacts of ``output_path``."""

    return output_path.with_name(f".{output_path.name}.tmp")


def _failure_sentinel_path(output_path: Path) -> Path:
    """Return the sentinel path recorded when a pipeline step fails."""

    return output_path.with_name(f"{output_path.name}.failed")


def _sidecar_paths(base: Path) -> dict[str, Path]:
    """Return known auxiliary artefacts generated from ``base``."""

    stemless = str(base.with_suffix(""))
    return {
        "meta": base.with_name(base.name + ".meta.yaml"),
        "failure": base.with_name(f"{base.stem}_failure_cases.csv"),
        "quality_json": base.with_suffix(".quality.json"),
        "quality_table": Path(f"{stemless}_quality_report_table.csv"),
        "corr_table": Path(f"{stemless}_data_correlation_report_table.csv"),
    }


def _run_step(
    step: PipelineStep,
    cfg: PipelineRunConfig,
    final_output: Path,
    working_output: Path,
) -> StepExecutionResult:
    """Execute ``step`` with ``cfg`` returning the resulting exit code."""

    input_path = step.required_input(cfg)
    if not input_path.exists():
        _LOGGER.error("step_input_missing", step=step.name, path=str(input_path))
        return StepExecutionResult(exit_code=1, executed=False)
    if cfg.skip_existing and final_output.exists() and not cfg.force:
        _LOGGER.info(
            "step_skipped_existing", step=step.name, path=str(final_output)
        )
        return StepExecutionResult(exit_code=0, executed=False)

    arguments = step.build_arguments(cfg, output_path=working_output)
    _LOGGER.debug("step_arguments", step=step.name, arguments=arguments)
    exit_code = step.main(arguments)
    return StepExecutionResult(exit_code=exit_code, executed=True)


def _finalize_step_success(
    final_output: Path, working_output: Path, sentinel_path: Path
) -> None:
    """Rename temporary outputs into place and clear failure sentinels."""

    final_sidecars = _sidecar_paths(final_output)
    working_sidecars = _sidecar_paths(working_output)

    if working_output.exists():
        if final_output.exists():
            final_output.unlink()
        working_output.replace(final_output)

    for name, working_path in working_sidecars.items():
        if not working_path.exists():
            continue
        final_path = final_sidecars[name]
        if final_path.exists():
            final_path.unlink()
        working_path.replace(final_path)
    if sentinel_path.exists():
        sentinel_path.unlink()


def _cleanup_failed_step(
    final_output: Path,
    working_output: Path,
    sentinel_path: Path,
    *,
    executed: bool,
) -> None:
    """Remove partial outputs and persist a failure sentinel."""

    final_sidecars = _sidecar_paths(final_output)
    working_sidecars = _sidecar_paths(working_output)

    candidates = [working_output]
    if executed:
        candidates.append(final_output)
    for candidate in candidates:
        if candidate.exists():
            try:
                candidate.unlink()
            except OSError as exc:  # pragma: no cover - defensive guard
                _LOGGER.warning(
                    "step_cleanup_failed",
                    path=str(candidate),
                    error=str(exc),
                )

    for name, working_path in working_sidecars.items():
        final_path = final_sidecars[name]
        if name == "failure":
            if working_path.exists():
                try:
                    if final_path.exists():
                        final_path.unlink()
                except OSError as exc:  # pragma: no cover - defensive guard
                    _LOGGER.warning(
                        "failure_sidecar_cleanup_failed",
                        path=str(final_path),
                        error=str(exc),
                    )
                try:
                    working_path.replace(final_path)
                except OSError as exc:  # pragma: no cover - defensive guard
                    _LOGGER.warning(
                        "failure_sidecar_promote_failed",
                        path=str(working_path),
                        error=str(exc),
                    )
            elif executed and final_path.exists():
                try:
                    final_path.unlink()
                except OSError as exc:  # pragma: no cover - defensive guard
                    _LOGGER.warning(
                        "failure_sidecar_cleanup_failed",
                        path=str(final_path),
                        error=str(exc),
                    )
            continue

        paths_to_remove = [working_path]
        if executed:
            paths_to_remove.append(final_path)
        for path in paths_to_remove:
            if path.exists():
                try:
                    path.unlink()
                except OSError as exc:  # pragma: no cover - defensive guard
                    _LOGGER.warning(
                        "step_cleanup_failed",
                        path=str(path),
                        error=str(exc),
                    )
    try:
        sentinel_path.touch()
    except OSError as exc:  # pragma: no cover - defensive guard
        _LOGGER.warning(
            "sentinel_write_failed",
            path=str(sentinel_path),
            error=str(exc),
        )


def run_pipeline(cfg: PipelineRunConfig) -> int:
    """Execute all configured steps and return the resulting exit status."""

    overall_status = 0
    _LOGGER.info("pipeline_start", stage="pipeline")
    for step in _PIPELINE_STEPS:
        _LOGGER.info("step_start", step=step.name)
        final_output = step.expected_output(cfg)
        working_output = _temporary_output_path(final_output)
        sentinel_path = _failure_sentinel_path(final_output)
        if working_output.exists():
            working_output.unlink()
        try:
            result = _run_step(step, cfg, final_output, working_output)
        except Exception as exc:  # pragma: no cover - defensive guard
            _LOGGER.exception("step_exception", step=step.name, error=str(exc))
            _cleanup_failed_step(
                final_output,
                working_output,
                sentinel_path,
                executed=True,
            )
            overall_status = 1
            break
        if result.exit_code != 0:
            _LOGGER.error(
                "step_failed", step=step.name, exit_code=result.exit_code
            )
            _cleanup_failed_step(
                final_output,
                working_output,
                sentinel_path,
                executed=result.executed,
            )
            overall_status = result.exit_code
            break
        _finalize_step_success(final_output, working_output, sentinel_path)
        _LOGGER.info("step_done", step=step.name)
    else:
        _LOGGER.info("workflow_complete")
    _LOGGER.info("pipeline_done", stage="pipeline", exit_code=overall_status)
    return overall_status


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry point for the orchestration script."""

    args = _parse_args(argv)
    try:
        logger = _configure_logging(args.log_level)
    except ValueError as exc:
        raise SystemExit(str(exc))

    global _LOGGER
    _LOGGER = logger
    try:
        cfg = _prepare_config(args)
    except (FileNotFoundError, OSError, ValueError) as exc:
        _LOGGER.error("configuration_error", error=str(exc))
        return 1
    status = run_pipeline(cfg)
    if status != 0:
        _LOGGER.error("workflow_failed", exit_code=status)
    else:
        _LOGGER.info("workflow_succeeded")
    return status


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

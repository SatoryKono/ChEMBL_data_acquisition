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

_PROJECT_ROOT = Path(__file__).resolve().parents[1]

if __package__ in {None, ""}:
    project_root_str = str(_PROJECT_ROOT)
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)

from library.utils.bootstrap import ensure_project_root


if __package__ in {None, ""}:
    ensure_project_root()

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
        output_csv = (
            output_path if output_path is not None else cfg.output_path(self.name)
        )
        args = ["--config", str(cfg.config_path), "--input", str(input_csv)]
        args.extend(["--output", str(output_csv)])
        args.extend(["--log-level", cfg.log_level])
        if cfg.limit is not None:
            args.extend(["--limit", str(cfg.limit)])
        if cfg.force:
            args.append("--force")
        if cfg.skip_existing:
            args.append("--skip-existing")
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
        help=(
            "Maximum number of identifiers processed by each pipeline; use 0 to "
            "skip execution"
        ),
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

    if args.limit is not None and args.limit < 0:
        raise ValueError("--limit must be zero or a positive integer")

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


def _coerce_exit_code(value: object) -> int:
    """Translate ``SystemExit.code`` values into integer exit codes."""

    if isinstance(value, int):
        return value
    if value is None:
        return 0
    if isinstance(value, str):
        try:
            return int(value)
        except ValueError:
            return 1
    return 1


@dataclass
class SidecarArtefact:
    """Describe auxiliary files associated with a pipeline output."""

    destination: Path
    final_path: Path | None = None
    working_path: Path | None = None


def _discover_sidecars(final_output: Path, working_output: Path) -> dict[Path, SidecarArtefact]:
    """Return all auxiliary files derived from ``final_output`` and ``working_output``."""

    final_dir = final_output.parent
    working_dir = working_output.parent
    sentinel_name = _failure_sentinel_path(final_output).name
    patterns = {
        final_output.name,
        final_output.with_suffix("").name,
        working_output.name,
        working_output.with_suffix("").name,
    }
    replacements = (
        (working_output.name, final_output.name),
        (working_output.with_suffix("").name, final_output.with_suffix("").name),
    )
    main_outputs = {final_output.resolve(), working_output.resolve()}

    def _normalise_relative(path: Path) -> Path:
        parts: list[str] = []
        for part in path.parts:
            normalised = part
            for old, new in replacements:
                if old == new or old not in normalised:
                    continue
                normalised = normalised.replace(old, new)
            parts.append(normalised)
        return Path(*parts)

    def _collect(base_dir: Path) -> dict[Path, Path]:
        collected: dict[Path, Path] = {}
        if not base_dir.exists():
            return collected
        for candidate in base_dir.rglob("*"):
            if not candidate.is_file():
                continue
            try:
                resolved = candidate.resolve()
            except OSError:  # pragma: no cover - filesystem race protection
                continue
            if resolved in main_outputs:
                continue
            name = candidate.name
            if name == sentinel_name:
                continue
            if not any(name.startswith(pattern) for pattern in patterns):
                continue
            rel_path = candidate.relative_to(base_dir)
            collected[rel_path] = candidate
        return collected

    sidecars: dict[Path, SidecarArtefact] = {}

    for rel_path, path in _collect(final_dir).items():
        canonical_rel = _normalise_relative(rel_path)
        destination = final_dir / canonical_rel
        entry = sidecars.get(canonical_rel)
        if entry is None:
            entry = SidecarArtefact(destination=destination)
            sidecars[canonical_rel] = entry
        entry.final_path = path

    for rel_path, path in _collect(working_dir).items():
        canonical_rel = _normalise_relative(rel_path)
        destination = final_dir / canonical_rel
        entry = sidecars.get(canonical_rel)
        if entry is None:
            entry = SidecarArtefact(destination=destination)
            sidecars[canonical_rel] = entry
        entry.working_path = path

    return sidecars


def _cleanup_empty_directories(path: Path, *, root: Path) -> None:
    """Remove empty directories upward from ``path`` until reaching ``root``."""

    try:
        root_resolved = root.resolve()
    except OSError:  # pragma: no cover - path vanished concurrently
        return
    current = path
    while True:
        try:
            current_resolved = current.resolve()
        except FileNotFoundError:
            parent = current.parent
            if parent == current:
                break
            current = parent
            continue
        except OSError:  # pragma: no cover - defensive guard
            break
        if current_resolved == root_resolved or current == root:
            break
        try:
            current.rmdir()
        except OSError:
            break
        parent = current.parent
        if parent == current:
            break
        current = parent


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
        _LOGGER.info("step_skipped_existing", step=step.name, path=str(final_output))
        return StepExecutionResult(exit_code=0, executed=False)

    arguments = step.build_arguments(cfg, output_path=working_output)
    _LOGGER.debug("step_arguments", step=step.name, arguments=arguments)
    try:
        exit_code = step.main(arguments)
    except SystemExit as exc:
        exit_code = _coerce_exit_code(exc.code)
        return StepExecutionResult(exit_code=exit_code, executed=True)
    except BaseException:
        raise
    return StepExecutionResult(exit_code=exit_code, executed=True)


def _finalize_step_success(
    final_output: Path, working_output: Path, sentinel_path: Path
) -> None:
    """Rename temporary outputs into place and clear failure sentinels."""

    sidecars = _discover_sidecars(final_output, working_output)
    working_dir = working_output.parent
    final_dir = final_output.parent

    if working_output.exists():
        if final_output.exists():
            final_output.unlink()
        working_output.replace(final_output)

    for sidecar in sidecars.values():
        working_path = sidecar.working_path
        if working_path is None or not working_path.exists():
            continue
        destination = sidecar.destination
        destination.parent.mkdir(parents=True, exist_ok=True)
        try:
            same_location = working_path.resolve() == destination.resolve()
        except OSError:
            same_location = False
        if same_location:
            final_path = sidecar.final_path
            if (
                final_path is not None
                and final_path.exists()
                and final_path != destination
            ):
                final_parent = final_path.parent
                final_path.unlink()
                _cleanup_empty_directories(final_parent, root=final_dir)
            continue
        if destination.exists():
            destination.unlink()
        original_parent = working_path.parent
        working_path.replace(destination)
        _cleanup_empty_directories(original_parent, root=working_dir)
        final_path = sidecar.final_path
        if (
            final_path is not None
            and final_path.exists()
            and final_path != destination
        ):
            final_parent = final_path.parent
            final_path.unlink()
            _cleanup_empty_directories(final_parent, root=final_dir)
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

    sidecars = _discover_sidecars(final_output, working_output)
    working_dir = working_output.parent
    final_dir = final_output.parent

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

    for sidecar in sidecars.values():
        destination = sidecar.destination
        is_failure = destination.name.endswith("_failure_cases.csv")
        working_path = sidecar.working_path
        final_path = sidecar.final_path

        if is_failure:
            if working_path is not None and working_path.exists():
                destination.parent.mkdir(parents=True, exist_ok=True)
                try:
                    if destination.exists():
                        destination.unlink()
                except OSError as exc:  # pragma: no cover - defensive guard
                    _LOGGER.warning(
                        "failure_sidecar_cleanup_failed",
                        path=str(destination),
                        error=str(exc),
                    )
                try:
                    original_parent = working_path.parent
                    working_path.replace(destination)
                    _cleanup_empty_directories(original_parent, root=working_dir)
                except OSError as exc:  # pragma: no cover - defensive guard
                    _LOGGER.warning(
                        "failure_sidecar_promote_failed",
                        path=str(working_path),
                        error=str(exc),
                    )
            elif executed and destination.exists():
                try:
                    destination.unlink()
                except OSError as exc:  # pragma: no cover - defensive guard
                    _LOGGER.warning(
                        "failure_sidecar_cleanup_failed",
                        path=str(destination),
                        error=str(exc),
                    )
                else:
                    _cleanup_empty_directories(destination.parent, root=final_dir)
            if (
                final_path is not None
                and final_path.exists()
                and final_path != destination
            ):
                try:
                    final_parent = final_path.parent
                    final_path.unlink()
                    _cleanup_empty_directories(final_parent, root=final_dir)
                except OSError as exc:  # pragma: no cover - defensive guard
                    _LOGGER.warning(
                        "failure_sidecar_cleanup_failed",
                        path=str(final_path),
                        error=str(exc),
                    )
            continue

        if working_path is not None and working_path.exists():
            try:
                original_parent = working_path.parent
                working_path.unlink()
                _cleanup_empty_directories(original_parent, root=working_dir)
            except OSError as exc:  # pragma: no cover - defensive guard
                _LOGGER.warning(
                    "step_cleanup_failed",
                    path=str(working_path),
                    error=str(exc),
                )
        if executed:
            removal_targets = {destination}
            if final_path is not None:
                removal_targets.add(final_path)
            for path in removal_targets:
                if path.exists():
                    try:
                        parent = path.parent
                        path.unlink()
                        _cleanup_empty_directories(parent, root=final_dir)
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
        except SystemExit as exc:  # pragma: no cover - defensive guard
            exit_code = _coerce_exit_code(exc.code)
            _LOGGER.error("step_system_exit", step=step.name, exit_code=exit_code)
            _cleanup_failed_step(
                final_output,
                working_output,
                sentinel_path,
                executed=True,
            )
            overall_status = exit_code
            break
        except BaseException as exc:  # pragma: no cover - defensive guard
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
            _LOGGER.error("step_failed", step=step.name, exit_code=result.exit_code)
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

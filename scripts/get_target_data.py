"""CLI entry-point for the reproducible target pipeline."""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING

# ruff: noqa: E402 - bootstrap adjusts import order for direct execution
if TYPE_CHECKING:
    from . import _bootstrap as _bootstrap_module
elif __package__ in {None, ""}:
    import _bootstrap as _bootstrap_module  # pragma: no cover - CLI fallback
else:  # pragma: no cover - executed when imported as a package module
    from . import _bootstrap as _bootstrap_module

bootstrap_cli = _bootstrap_module.bootstrap_cli
bootstrap_cli(__package__, __file__)
del bootstrap_cli
del _bootstrap_module

from library.io.config_loader import load_config
from library.io.metadata_writer import save_metadata
from library.io.output_writer import save_standard_outputs
from library.postprocessing.target import steps as target_steps
from library.utils.logging import get_logger

TABLE_NAME = "target"


@dataclass(frozen=True)
class PipelineArgs:
    """Normalised CLI arguments for the target pipeline."""

    command: str
    limit: int
    date_tag: str
    output_dir: Path
    config: Path | None
    log_level: str
    input_dir: Path


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("limit must be positive")
    return parsed


def _log_level(value: str) -> str:
    """Return a normalised logging level name."""

    levels = {"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"}
    candidate = value.strip().upper()
    if candidate not in levels:
        options = ", ".join(sorted(levels))
        raise argparse.ArgumentTypeError(f"log level must be one of: {options}")
    return candidate


def _configure_logging(level_name: str) -> None:
    """Apply the requested log level to the root logger."""

    level = getattr(logging, level_name, logging.INFO)
    logging.getLogger().setLevel(level)


def _date_tag(value: str) -> str:
    """Validate that ``value`` is formatted as ``YYYYMMDD``."""

    try:
        datetime.strptime(value, "%Y%m%d")
    except ValueError as exc:  # noqa: TRY003 - convert into argparse error
        raise argparse.ArgumentTypeError("date tag must be formatted as YYYYMMDD") from exc
    return value


def _log_level(value: str) -> str:
    """Ensure ``value`` resolves to a valid logging level name or number."""

    candidate = value.strip()
    if not candidate:
        raise argparse.ArgumentTypeError("log level must not be empty")
    try:
        _resolve_log_level(candidate)
    except ValueError as exc:  # noqa: TRY003 - convert into argparse error
        raise argparse.ArgumentTypeError(str(exc)) from exc
    return candidate


def _resolve_log_level(value: str) -> int:
    """Return the numeric logging level represented by ``value``."""

    normalised = value.strip()
    if normalised.isdigit():
        numeric = int(normalised)
        if numeric < 0:
            raise ValueError("log level must be non-negative")
        return numeric
    resolved = logging.getLevelName(normalised.upper())
    if isinstance(resolved, int):
        return resolved
    raise ValueError(f"unknown log level: {value}")


def _configure_logging(value: str) -> None:
    """Configure the root logger using ``value`` parsed from CLI arguments."""

    logging.getLogger().setLevel(_resolve_log_level(value))


def parse_args(argv: Sequence[str] | None = None) -> tuple[argparse.Namespace, PipelineArgs]:
    """Parse CLI arguments returning both the raw namespace and dataclass view."""

    parser = argparse.ArgumentParser(description="Fetch and enrich ChEMBL target data")
    parser.add_argument(
        "command",
        nargs="?",
        default="all",
        help="Execution mode compatible with legacy target pipelines",
    )
    parser.add_argument("--limit", type=_positive_int, default=1000, help="Number of targets to fetch")
    parser.add_argument(
        "--date-tag",
        type=_date_tag,
        default=datetime.now(UTC).strftime("%Y%m%d"),
        help="Execution date token in YYYYMMDD format",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/output"),
        help="Directory for generated artefacts",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("data/input"),
        help="Directory containing cached inputs for compatibility",
    )
    parser.add_argument(
        "--log-level",
        type=_log_level,
        default="INFO",
        help="Logging level (numeric or name)",
    )
    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to configuration file (defaults to config/config.yaml)",
    )
    parser.add_argument(
        "--log-level",
        type=_log_level,
        default="INFO",
        help="Root logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("data/input"),
        help="Directory containing auxiliary input artefacts",
    )
    namespace = parser.parse_args(argv)
    config_path = Path(namespace.config) if namespace.config is not None else None
    pipeline_args = PipelineArgs(
        command=namespace.command,
        limit=namespace.limit,
        date_tag=namespace.date_tag,
        output_dir=Path(namespace.output_dir),
        config=config_path,
        log_level=namespace.log_level,
        input_dir=Path(namespace.input_dir),
    )
    return namespace, pipeline_args


def main(argv: Sequence[str] | None = None) -> int:
    """Execute the target pipeline."""

    namespace, args = parse_args(argv)
    _configure_logging(args.log_level)
    logger = get_logger(__name__)

    if args.command != "all":
        logger.error("target_unknown_command", command=args.command)
        return 2

    try:
        config = load_config(args.config) if args.config is not None else load_config(None)
    except Exception as exc:  # noqa: BLE001 - surface configuration failures
        logger.error("target_config_error", error=str(exc))
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.input_dir.mkdir(parents=True, exist_ok=True)

    try:
        dataset = target_steps.fetch_normalize_target(args.limit)
        target_data = target_steps.generate_target_reports(dataset)
        artifacts = save_standard_outputs(
            target_data.dataset,
            target_data.correlation_report,
            target_data.quality_report,
            TABLE_NAME,
            args.date_tag,
            output_dir=args.output_dir,
        )
        artifact_paths = [artifacts.dataset, artifacts.correlation_report, artifacts.quality_report]
        save_metadata(
            TABLE_NAME,
            args.date_tag,
            namespace,
            qc_summary=target_data.qc_summary,
            output_dir=args.output_dir,
            artifacts=artifact_paths,
            sources=[config.path.name],
        )
    except Exception as exc:  # noqa: BLE001 - top-level failure boundary
        logger.exception("target_pipeline_failed", error=str(exc))
        return 1

    logger.info("target_pipeline_success", rows=target_data.dataset.shape[0])
    return 0


if __name__ == "__main__":
    sys.exit(main())

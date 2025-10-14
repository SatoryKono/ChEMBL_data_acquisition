"""CLI entry-point for the reproducible target pipeline."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

if __package__ in {None, ""}:
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

from library.io.config_loader import load_config
from library.io.metadata_writer import save_metadata
from library.io.output_writer import save_standard_outputs
from library.postprocessing.target.steps import (
    fetch_normalize_target,
    generate_target_reports,
)
from library.utils.logging import get_logger

TABLE_NAME = "target"


@dataclass(frozen=True)
class PipelineArgs:
    """Normalised CLI arguments for the target pipeline."""

    limit: int
    date_tag: str
    output_dir: Path
    config: Path | None


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("limit must be positive")
    return parsed


def _date_tag(value: str) -> str:
    try:
        datetime.strptime(value, "%Y%m%d")
    except ValueError as exc:  # noqa: TRY003 - convert into argparse error
        raise argparse.ArgumentTypeError("date tag must be formatted as YYYYMMDD") from exc
    return value


def parse_args(argv: Sequence[str] | None = None) -> tuple[argparse.Namespace, PipelineArgs]:
    """Parse CLI arguments returning both the raw namespace and dataclass view."""

    parser = argparse.ArgumentParser(description="Fetch and enrich ChEMBL target data")
    parser.add_argument("--limit", type=_positive_int, default=1000, help="Number of targets to fetch")
    parser.add_argument(
        "--date-tag",
        type=_date_tag,
        default=datetime.utcnow().strftime("%Y%m%d"),
        help="Execution date token in YYYYMMDD format",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/output"),
        help="Directory for generated artefacts",
    )
    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to configuration file (defaults to config/config.yaml)",
    )
    namespace = parser.parse_args(argv)
    config_path = Path(namespace.config) if namespace.config is not None else None
    pipeline_args = PipelineArgs(
        limit=namespace.limit,
        date_tag=namespace.date_tag,
        output_dir=Path(namespace.output_dir),
        config=config_path,
    )
    return namespace, pipeline_args


def main(argv: Sequence[str] | None = None) -> int:
    """Execute the target pipeline."""

    namespace, args = parse_args(argv)
    logger = get_logger(__name__)

    try:
        config = load_config(args.config) if args.config is not None else load_config(None)
    except Exception as exc:  # noqa: BLE001 - surface configuration failures
        logger.error("target_config_error", error=str(exc))
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)

    try:
        dataset = fetch_normalize_target(args.limit)
        target_data = generate_target_reports(dataset)
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

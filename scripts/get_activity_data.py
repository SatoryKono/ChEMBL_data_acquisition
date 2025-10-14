"""CLI entry point producing normalized ChEMBL activity exports."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Sequence

from library.io.config_loader import Config, load_config
from library.io.output_writer import save_standard_outputs
from library.postprocessing.activity.steps import (
    ACTIVITY_COLUMN_ORDER,
    ACTIVITY_SORT_COLUMNS,
    run_activity_pipeline,
)
from library.utils.logging import Logger, get_logger

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "data" / "output"
DEFAULT_CONFIG_PATH = Path("config/config.yaml")


@dataclass(frozen=True)
class PipelineArgs:
    """Container describing the normalized activity pipeline invocation."""

    limit: int | None
    date_tag: str
    output_dir: Path
    config_path: Path


def _default_date_tag() -> str:
    return datetime.now(tz=UTC).strftime("%Y%m%d")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=None, help="Maximum number of rows to fetch")
    parser.add_argument(
        "--date-tag",
        dest="date_tag",
        default=None,
        help="Date token used in output artefact names (YYYYMMDD)",
    )
    parser.add_argument(
        "--output-dir",
        dest="output_dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Destination directory for generated artefacts",
    )
    parser.add_argument(
        "--config",
        dest="config_path",
        default=str(DEFAULT_CONFIG_PATH),
        help="Configuration file controlling API parameters",
    )
    return parser


def _normalize_output_dir(value: str | Path) -> Path:
    candidate = Path(value)
    if candidate.is_absolute():
        return candidate
    return (ROOT / candidate).resolve()


def _normalize_config_path(value: str | Path) -> Path:
    candidate = Path(value)
    if candidate.is_absolute():
        return candidate
    return ROOT / candidate


def parse_args(argv: Sequence[str] | None = None) -> PipelineArgs:
    parser = build_parser()
    namespace = parser.parse_args(list(argv) if argv is not None else None)
    limit = namespace.limit
    if limit is not None and limit < 0:
        parser.error("--limit must be non-negative")
    date_tag = namespace.date_tag or _default_date_tag()
    output_dir = _normalize_output_dir(namespace.output_dir)
    config_path = _normalize_config_path(namespace.config_path)
    return PipelineArgs(limit=limit, date_tag=date_tag, output_dir=output_dir, config_path=config_path)


def _resolve_chembl_base(cfg: Config) -> str | None:
    data = cfg.to_dict()
    sources = data.get("sources")
    if not isinstance(sources, dict):
        return None
    chembl_cfg = sources.get("chembl")
    if not isinstance(chembl_cfg, dict):
        return None
    api_cfg = chembl_cfg.get("api")
    if not isinstance(api_cfg, dict):
        return None
    base = api_cfg.get("chembl_base")
    if isinstance(base, str) and base.strip():
        return base.strip()
    return None


def _prepare_logger(pipeline_args: PipelineArgs) -> Logger:
    logger = get_logger(__name__)
    return logger.bind(
        limit=pipeline_args.limit,
        output_dir=str(pipeline_args.output_dir),
        date_tag=pipeline_args.date_tag,
    )


def run_pipeline(pipeline_args: PipelineArgs, *, logger: Logger | None = None) -> int:
    log = logger or _prepare_logger(pipeline_args)
    log.info("activity_pipeline_start")
    try:
        cfg = load_config(pipeline_args.config_path)
    except Exception as exc:  # pragma: no cover - configuration errors are fatal
        log.error(
            "activity_config_load_failed",
            error=str(exc),
            config=str(pipeline_args.config_path),
        )
        return 1

    base_url = _resolve_chembl_base(cfg)
    dataset, correlation, quality = run_activity_pipeline(
        limit=pipeline_args.limit,
        base_url=base_url,
        logger=log.bind(stage="activity_pipeline"),
    )

    artifacts = save_standard_outputs(
        dataset,
        correlation,
        quality,
        table_name="activity",
        date_tag=pipeline_args.date_tag,
        key_columns=ACTIVITY_SORT_COLUMNS,
        column_order=ACTIVITY_COLUMN_ORDER,
        output_dir=pipeline_args.output_dir,
    )

    log.info(
        "activity_pipeline_complete",
        rows=int(dataset.shape[0]),
        dataset=str(artifacts.dataset),
        quality=str(artifacts.quality_report),
        correlation=str(artifacts.correlation_report),
    )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    pipeline_args = parse_args(argv)
    return run_pipeline(pipeline_args)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

"""Helpers for persisting pipeline metadata sidecars."""

from __future__ import annotations

import argparse
import logging
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml


logger = logging.getLogger(__name__)


def _serialise_value(value: Any) -> Any:
    """Return a YAML-safe representation of ``value``."""

    if isinstance(value, Path):
        return str(value)
    if isinstance(value, argparse.Namespace):
        return {key: _serialise_value(val) for key, val in vars(value).items()}
    if isinstance(value, Mapping):
        return {str(key): _serialise_value(val) for key, val in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_serialise_value(item) for item in value]
    if isinstance(value, datetime):
        if value.tzinfo is None:
            return value.replace(tzinfo=timezone.utc).isoformat().replace("+00:00", "Z")
        return value.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")
    return value


def _serialise_parameters(args: argparse.Namespace | Mapping[str, Any] | None) -> dict[str, Any]:
    """Convert CLI arguments to a mapping suitable for YAML serialisation."""

    if args is None:
        return {}

    if isinstance(args, Mapping):
        source = dict(args)
    else:
        source = vars(args).copy()

    return {key: _serialise_value(value) for key, value in source.items()}


def _resolve_outputs(
    table_name: str,
    date_tag: str,
    artifacts: Sequence[Path] | None = None,
) -> list[str]:
    """Return the list of output artefact names recorded in metadata."""

    if artifacts:
        return [str(path.name) for path in artifacts]

    stem = f"output.{table_name}_{date_tag}"
    return [
        f"{stem}.csv",
        f"{stem}_quality_report_table.csv",
        f"{stem}_data_correlation_report_table.csv",
    ]


def save_metadata(
    table_name: str,
    date_tag: str,
    args: argparse.Namespace | Mapping[str, Any] | None,
    *,
    qc_summary: Mapping[str, Any] | None = None,
    output_dir: Path | str | None = None,
    artifacts: Sequence[Path] | None = None,
    sources: Sequence[str] | None = None,
) -> Path:
    """Persist pipeline metadata for ``table_name`` outputs."""

    resolved_output_dir = Path(output_dir) if output_dir is not None else Path("data/output")
    resolved_output_dir.mkdir(parents=True, exist_ok=True)

    meta_path = resolved_output_dir / f"output.{table_name}_{date_tag}.meta.yaml"

    parameters = _serialise_parameters(args)
    outputs = _resolve_outputs(table_name, date_tag, artifacts=artifacts)

    meta: dict[str, Any] = {
        "table": table_name,
        "generated_at": datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "pipeline_version": "2.1",
        "parameters": parameters,
        "sources": list(sources or []),
        "outputs": outputs,
        "qc_summary": dict(qc_summary or {}),
    }

    with meta_path.open("w", encoding="utf-8") as buffer:
        yaml.safe_dump(meta, buffer, sort_keys=False, allow_unicode=True)

    logger.info("[META] Метаданные сохранены: %s", meta_path)
    return meta_path


__all__ = ["save_metadata"]

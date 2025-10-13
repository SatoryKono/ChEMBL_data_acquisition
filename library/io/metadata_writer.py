"""Helpers for persisting pipeline metadata sidecars."""

from __future__ import annotations

import argparse
import logging
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from dataclasses import asdict, is_dataclass

from library.common.metadata_writer import write_meta_yaml
from library.common.run_context import RunContext, get_current


logger = logging.getLogger(__name__)


def _serialise_value(value: Any) -> Any:
    """Return a YAML-safe representation of ``value``."""

    if is_dataclass(value):
        return _serialise_value(asdict(value))
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
    if callable(value):
        module = getattr(value, "__module__", None)
        qualname = getattr(value, "__qualname__", None)
        if module and qualname:
            return f"{module}.{qualname}"
        return repr(value)
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


def _resolve_generated_at(run_context: RunContext | None) -> str:
    """Return the timestamp to be stored in metadata."""

    if run_context is not None:
        generated_at = getattr(run_context, "generated_at", None)
        if generated_at:
            return str(generated_at)

    return (
        datetime.utcnow()
        .replace(tzinfo=timezone.utc)
        .isoformat(timespec="seconds")
        .replace("+00:00", "Z")
    )


def save_metadata(
    table_name: str,
    date_tag: str,
    args: argparse.Namespace | Mapping[str, Any] | None,
    *,
    qc_summary: Mapping[str, Any] | None = None,
    output_dir: Path | str | None = None,
    artifacts: Sequence[Path] | None = None,
    sources: Sequence[str] | None = None,
    run_context: RunContext | None = None,
    stats_extra: Mapping[str, Any] | None = None,
) -> Path:
    """Persist pipeline metadata for ``table_name`` outputs."""

    resolved_output_dir = Path(output_dir) if output_dir is not None else Path("data/output")

    parameters = _serialise_parameters(args)
    outputs = _resolve_outputs(table_name, date_tag, artifacts=artifacts)

    generated_at = datetime.now(timezone.utc).isoformat(timespec="seconds").replace(
        "+00:00",
        "Z",
    )

    csv_stub_path = resolved_output_dir / f"output.{table_name}_{date_tag}"

    meta_path = write_meta_yaml(
        csv_path=csv_stub_path,
        command="",
        generated_at=generated_at,
        stats=dict(stats_extra) if stats_extra else None,
        extra_metadata={
            "table": table_name,
            "parameters": parameters,
            "sources": list(sources or []),
            "outputs": outputs,
            "qc_summary": dict(qc_summary or {}),
        },
    )

    logger.info("[META] Метаданные сохранены: %s", meta_path)
    return meta_path


__all__ = ["save_metadata"]

"""Helpers for producing standardised pipeline reporting artefacts."""

from __future__ import annotations

import json
from collections.abc import Callable, Mapping, MutableMapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import yaml

from ..common.metadata import Stats, file_sha256, write_meta_yaml
from ..qa.reporting import TableQualityHook
from ..qa.table_quality import TableQualityProfiler, analyze_table_quality


class RunManifestError(RuntimeError):
    """Base exception for run manifest reporting failures."""


class QualityReportError(RunManifestError):
    """Raised when the quality report cannot be serialised to disk."""

    def __init__(self, message: str, path: Path | None = None) -> None:
        super().__init__(message)
        self.path = Path(path) if path is not None else None


class QualityAnalysisError(RunManifestError):
    """Raised when table quality analysis fails."""


@dataclass(frozen=True, slots=True)
class PipelineOutputReport:
    """Container capturing output artefact metadata for a pipeline step."""

    csv_path: Path
    stats: Stats
    meta_path: Path
    meta_sha256: str
    quality_path: Path | None = None
    quality_sha256: str | None = None

    def as_manifest_payload(self) -> dict[str, Any]:
        """Return a serialisable payload summarising the artefact."""

        stats_payload: dict[str, Any] = {}
        for key, value in self.stats.items():
            if key in {"rows_total", "rows_kept", "rows_dropped"}:
                stats_payload[key] = int(value) if value is not None else None
            else:
                stats_payload[key] = value
        payload: dict[str, Any] = {
            "output_path": str(self.csv_path),
            "output_sha256": self.stats.get("output_sha256"),
            "stats": stats_payload,
            "meta": {
                "path": str(self.meta_path),
                "checksum_sha256": self.meta_sha256,
            },
        }
        if self.quality_path is not None:
            payload["quality"] = {
                "path": str(self.quality_path),
                "checksum_sha256": self.quality_sha256,
            }
        return payload


def _build_stats(
    csv_path: Path,
    *,
    rows_total: int,
    rows_kept: int,
    extra: Mapping[str, Any] | None = None,
) -> Stats:
    total = int(rows_total)
    kept = int(rows_kept)
    dropped = max(total - kept, 0)
    stats: dict[str, Any] = {
        "rows_total": total,
        "rows_kept": kept,
        "rows_dropped": dropped,
        "output_sha256": file_sha256(csv_path),
    }
    if extra:
        stats.update(extra)
    return cast(Stats, stats)


def _normalise_quality_summary(
    summary: Any,
    builder: Callable[[Any], Mapping[str, Any]] | None,
) -> Mapping[str, Any]:
    if builder is not None:
        report = builder(summary)
        if not isinstance(report, Mapping):
            raise TypeError("quality report builder must return a mapping")
        return dict(report)
    if summary is None:
        raise TypeError("quality_summary is required when no builder is provided")
    if isinstance(summary, Mapping):
        return dict(summary)
    build = getattr(summary, "build", None)
    if callable(build):
        built = build()
        if isinstance(built, Mapping):
            return dict(built)
    raise TypeError("quality_summary must be a mapping or provide a builder")


def _cfg_value(cfg: Any, key: str, default: Any = None) -> Any:
    if cfg is None:
        return default
    if isinstance(cfg, Mapping):
        return cfg.get(key, default)
    return getattr(cfg, key, default)


def finalise_csv_output(
    *,
    csv_path: Path,
    rows_total: int,
    rows_kept: int,
    command: str,
    config: Mapping[str, Any],
    inputs: Mapping[str, Any],
    schema: str,
    stats_extra: Mapping[str, Any] | None = None,
    invocation: Sequence[str] | None = None,
    extra_metadata: Mapping[str, Any] | None = None,
    quality_summary: Any | None = None,
    quality_builder: Callable[[Any], Mapping[str, Any]] | None = None,
    quality_path: Path | None = None,
    quality_profiler: TableQualityProfiler | None = None,
    quality_config: Any | None = None,
    quality_table_name: str | None = None,
    quality_destination: Path | None = None,
    quality_hook: TableQualityHook | None = None,
) -> PipelineOutputReport:
    """Write metadata and optional quality artefacts for ``csv_path``."""

    path = Path(csv_path)
    stats = _build_stats(
        path, rows_total=rows_total, rows_kept=rows_kept, extra=stats_extra
    )
    meta_path = write_meta_yaml(
        csv_path=path,
        command=command,
        config=config,
        inputs=inputs,
        stats=stats,
        schema=schema,
        invocation=invocation,
        extra_metadata=extra_metadata,
    )
    meta_sha = file_sha256(meta_path)

    resolved_quality_path: Path | None = None
    quality_sha: str | None = None
    if quality_summary is not None or quality_builder is not None:
        report_data = _normalise_quality_summary(quality_summary, quality_builder)
        destination = (
            Path(quality_path)
            if quality_path is not None
            else path.with_suffix(".quality.json")
        )
        try:
            destination.write_text(
                json.dumps(report_data, ensure_ascii=False, indent=2),
                encoding="utf-8",
            )
        except (OSError, TypeError, ValueError) as exc:  # pragma: no cover - defensive
            raise QualityReportError(str(exc), destination) from exc
        resolved_quality_path = destination
        quality_sha = file_sha256(destination)

    quality_subject = quality_profiler if quality_profiler is not None else path
    if quality_hook is not None:
        try:
            quality_hook(quality_subject)
        except Exception as exc:  # pragma: no cover - propagated to caller
            raise QualityAnalysisError(str(exc)) from exc
    elif quality_profiler is not None and bool(
        _cfg_value(quality_config, "enable", False)
    ):
        table_name = quality_table_name or path.with_suffix("").name
        destination_dir = quality_destination or path.parent
        try:
            analyze_table_quality(
                quality_profiler,
                table_name=table_name,
                destination_dir=destination_dir,
                sample_rows=_cfg_value(quality_config, "sample_rows"),
                include_columns=_cfg_value(quality_config, "include_columns"),
                exclude_columns=_cfg_value(quality_config, "exclude_columns"),
            )
        except Exception as exc:  # pragma: no cover - propagated to caller
            raise QualityAnalysisError(str(exc)) from exc

    return PipelineOutputReport(
        csv_path=path,
        stats=stats,
        meta_path=meta_path,
        meta_sha256=meta_sha,
        quality_path=resolved_quality_path,
        quality_sha256=quality_sha,
    )


def merge_run_output(
    entry: MutableMapping[str, Any], report: PipelineOutputReport
) -> None:
    """Merge ``report`` statistics into a manifest ``entry`` in place."""

    stats: dict[str, Any] = {}
    for key, value in report.stats.items():
        if key in {"rows_total", "rows_kept", "rows_dropped"}:
            try:
                stats[key] = int(value)
            except (TypeError, ValueError):
                stats[key] = value
        else:
            stats[key] = value
    entry["stats"] = stats

    output = entry.setdefault("output", {})
    output.setdefault("path", str(report.csv_path))
    output["checksum_sha256"] = report.stats.get("output_sha256")
    output["meta_path"] = str(report.meta_path)
    output["meta_sha256"] = report.meta_sha256
    if report.quality_path is not None:
        output["quality_path"] = str(report.quality_path)
        output["quality_sha256"] = report.quality_sha256


def load_output_report(csv_path: Path) -> PipelineOutputReport | None:
    """Return a :class:`PipelineOutputReport` derived from metadata files."""

    path = Path(csv_path)
    meta_path = path.with_name(path.name + ".meta.yaml")
    if not meta_path.exists():
        return None
    try:
        with meta_path.open("r", encoding="utf-8") as fh:
            meta_data = yaml.safe_load(fh) or {}
    except (OSError, yaml.YAMLError):  # pragma: no cover - defensive parsing
        return None
    stats_raw = meta_data.get("stats") if isinstance(meta_data, Mapping) else None
    if not isinstance(stats_raw, Mapping):
        return None
    stats: Stats = {
        "rows_total": int(stats_raw.get("rows_total", 0)),
        "rows_kept": int(stats_raw.get("rows_kept", 0)),
        "rows_dropped": int(stats_raw.get("rows_dropped", 0)),
        "output_sha256": str(stats_raw.get("output_sha256", "")),
    }
    for key, value in stats_raw.items():
        if key not in stats:
            stats[key] = value
    if not stats.get("output_sha256") and path.exists():
        stats["output_sha256"] = file_sha256(path)
    meta_sha = file_sha256(meta_path)

    quality_path = path.with_suffix(".quality.json")
    quality_sha: str | None = None
    if quality_path.exists():
        quality_sha = file_sha256(quality_path)
    return PipelineOutputReport(
        csv_path=path,
        stats=stats,
        meta_path=meta_path,
        meta_sha256=meta_sha,
        quality_path=quality_path if quality_path.exists() else None,
        quality_sha256=quality_sha,
    )


__all__ = [
    "PipelineOutputReport",
    "QualityAnalysisError",
    "QualityReportError",
    "RunManifestError",
    "finalise_csv_output",
    "load_output_report",
    "merge_run_output",
]

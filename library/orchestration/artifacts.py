"""Utilities for managing workflow artefacts and diagnostics."""

from __future__ import annotations

import logging
import time
from collections import deque
from dataclasses import dataclass
from fnmatch import fnmatch
from pathlib import Path
from typing import Any, Iterable, Sequence

from library.io.paths import derive_output_labels
from library.postprocess.common import (
    SUPPORTED_TABLES as POSTPROCESS_SUPPORTED_TABLES,
)
from library.postprocess.common import (
    PostprocessResult,
    PostprocessingPipelineConfig,
    get_csv_runtime_config as get_postprocess_csv_config,
    get_pipeline_config as load_postprocess_pipeline_config,
    run_postprocessing_pipeline,
)
from library.postprocessing.activities import ACTIVITY_SCHEMA, run_activity_pipeline
from library.postprocessing.activities import validate_activities as validate_activity_reports
from library.postprocessing.assays import ASSAY_SCHEMA, run_assay_pipeline
from library.postprocessing.assays import validate_assays as validate_assay_reports
from library.postprocessing.documents import DOCUMENT_SCHEMA, run_document_pipeline
from library.postprocessing.documents import validate_documents as validate_document_reports
from library.postprocessing.targets import TARGET_SCHEMA, run_target_pipeline
from library.postprocessing.targets import validate_targets as validate_target_reports
from library.postprocessing.testitem import (
    TESTITEM_SCHEMA,
    run_testitem_pipeline,
    validate_testitems,
)

_UNLINK_MAX_ATTEMPTS = 5
_UNLINK_RETRY_SLEEP_SECONDS = 0.2
_WINDOWS_SHARING_VIOLATION = 32


@dataclass
class SidecarArtefact:
    """Describe auxiliary files associated with a pipeline output."""

    destination: Path
    final_path: Path | None = None
    working_path: Path | None = None


def remove_path(path: Path, *, logger: logging.Logger | None = None) -> None:
    """Remove ``path`` handling transient Windows sharing violations."""

    if not path.exists():
        return
    for attempt in range(1, _UNLINK_MAX_ATTEMPTS + 1):
        try:
            path.unlink()
            return
        except FileNotFoundError:
            return
        except PermissionError as exc:
            if attempt == _UNLINK_MAX_ATTEMPTS:
                raise
            if logger is not None:
                logger.debug(
                    "unlink_retry_permission",
                    path=str(path),
                    attempt=attempt,
                    error=str(exc),
                    exc_info=exc,
                )
            time.sleep(_UNLINK_RETRY_SLEEP_SECONDS)
            continue
        except OSError as exc:
            winerror = getattr(exc, "winerror", None)
            if (
                winerror == _WINDOWS_SHARING_VIOLATION
                and attempt < _UNLINK_MAX_ATTEMPTS
            ):
                if logger is not None:
                    logger.debug(
                        "unlink_retry_sharing_violation",
                        path=str(path),
                        attempt=attempt,
                        error=str(exc),
                        exc_info=exc,
                    )
                time.sleep(_UNLINK_RETRY_SLEEP_SECONDS)
                continue
            raise


def cleanup_empty_directories(path: Path, *, root: Path) -> None:
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


def discover_sidecars(
    final_output: Path,
    working_output: Path,
    *,
    max_depth: int | None = None,
    include_patterns: Sequence[str] | None = None,
    sentinel_path: Path | None = None,
) -> dict[Path, SidecarArtefact]:
    """Return all auxiliary files derived from ``final_output`` and ``working_output``."""

    final_dir = final_output.parent
    working_dir = working_output.parent
    if sentinel_path is None:
        sentinel_path = final_output.with_name(f"{final_output.name}.failed")
    sentinel_name = sentinel_path.name
    patterns = {
        final_output.name,
        final_output.with_suffix("").name,
        working_output.name,
        working_output.with_suffix("").name,
    }
    prefix_candidates: set[str] = set()
    for candidate in (final_output, working_output):
        current = candidate
        while True:
            prefix_candidates.add(current.name)
            if not current.suffix:
                break
            current = current.with_suffix("")
    replacements = (
        (working_output.name, final_output.name),
        (working_output.with_suffix("").name, final_output.with_suffix("").name),
    )
    main_outputs = {final_output.resolve(), working_output.resolve()}

    def _matches_explicit_patterns(rel_path: Path, name: str) -> bool:
        if include_patterns is None:
            return False
        rel_value = rel_path.as_posix()
        return any(
            fnmatch(rel_value, pattern) or fnmatch(name, pattern)
            for pattern in include_patterns
        )

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
        queue: deque[tuple[Path, int, bool]] = deque([(base_dir, 0, False)])
        while queue:
            current, depth, within_branch = queue.popleft()
            try:
                entries = list(current.iterdir())
            except OSError:  # pragma: no cover - filesystem race protection
                continue
            for entry in entries:
                rel_path = entry.relative_to(base_dir)
                name = entry.name
                if entry.is_dir():
                    next_depth = depth + 1
                    if max_depth is not None and next_depth > max_depth:
                        continue
                    matches_prefix = any(
                        name.startswith(prefix) for prefix in prefix_candidates
                    )
                    matches_pattern = _matches_explicit_patterns(rel_path, name)
                    next_within = within_branch or matches_prefix or matches_pattern
                    if next_within:
                        queue.append((entry, next_depth, next_within))
                    continue
                matches_default = any(name.startswith(pattern) for pattern in patterns)
                matches_explicit = _matches_explicit_patterns(rel_path, name)
                if not matches_default and not matches_explicit:
                    continue
                if name == sentinel_name:
                    continue
                try:
                    resolved = entry.resolve()
                except OSError:  # pragma: no cover - filesystem race protection
                    continue
                if resolved in main_outputs:
                    continue
                collected[rel_path] = entry
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


def is_diagnostic_sidecar(path: Path) -> bool:
    """Return ``True`` when ``path`` points to a diagnostic-only artefact."""

    name = path.name.lower()
    if name.startswith("output.") and ".csv_" in name:
        return True
    if name.endswith(".meta.yaml"):
        return True
    if name.endswith(".postprocess.report.json"):
        return True
    if name.endswith(".quality.json"):
        return True
    if name.startswith("output_postprocessed.") and name.endswith(".csv"):
        return True

    diagnostic_suffixes = (
        "_failure_cases.csv",
        "_chembl.csv",
        "_pubmed.csv",
        "_uniprot.csv",
        "_iuphar.csv",
    )
    if name.endswith(diagnostic_suffixes):
        return True

    if name.endswith(".csv") and (
        "_normalized" in name or "_raw" in name or name.endswith("_raw")
    ):
        return True

    raw_suffixes = (
        ".raw.csv",
        ".raw.parquet",
        ".raw.json",
        ".raw.jsonl",
        "_raw.parquet",
        "_raw.json",
        "_raw.jsonl",
    )
    if any(name.endswith(suffix) for suffix in raw_suffixes):
        return True

    return False


def cleanup_diagnostics(
    sidecars: Iterable[SidecarArtefact],
    *,
    diagnostics_enabled: bool,
    working_root: Path,
    final_root: Path,
    canonical_prefix: str | None,
    canonical_allowed: frozenset[str],
    logger: logging.Logger,
) -> None:
    """Remove diagnostic artefacts when diagnostics are disabled."""

    if diagnostics_enabled:
        return

    for sidecar in sidecars:
        working_path = sidecar.working_path
        destination = sidecar.destination
        extra_standard_output = False
        if canonical_prefix is not None:
            candidate_name = destination.name
            if (
                candidate_name.startswith(canonical_prefix)
                and candidate_name not in canonical_allowed
            ):
                extra_standard_output = True
        if not is_diagnostic_sidecar(destination) and not extra_standard_output:
            continue
        for candidate in (working_path, sidecar.final_path, destination):
            if candidate is None or not candidate.exists():
                continue
            try:
                remove_path(candidate, logger=logger)
            except OSError as exc:  # pragma: no cover - defensive guard
                logger.warning(
                    "diagnostic_sidecar_cleanup_failed",
                    path=str(candidate),
                    error=str(exc),
                )
                continue
            parent = candidate.parent
            root = working_root if parent.is_relative_to(working_root) else final_root
            cleanup_empty_directories(parent, root=root)


def finalize_outputs(
    final_output: Path,
    working_output: Path,
    sentinel_path: Path,
    *,
    diagnostics_enabled: bool,
    logger: logging.Logger,
    include_patterns: Sequence[str] | None = None,
) -> None:
    """Rename temporary outputs into place and clear failure sentinels."""

    sidecars = discover_sidecars(
        final_output,
        working_output,
        include_patterns=include_patterns,
    )
    working_dir = working_output.parent
    final_dir = final_output.parent

    try:
        table_name, date_tag = derive_output_labels(
            final_output.name,
            default_table=final_output.stem or "dataset",
        )
    except Exception:  # pragma: no cover - defensive guard
        canonical_prefix: str | None = None
        canonical_allowed: frozenset[str] = frozenset()
    else:
        canonical_prefix = f"output.{table_name}_{date_tag}"
        canonical_allowed = frozenset(
            {
                f"{canonical_prefix}.csv",
                f"{canonical_prefix}_data_correlation_report_table.csv",
                f"{canonical_prefix}_quality_report_table.csv",
            }
        )

    if working_output.exists():
        if final_output.exists():
            remove_path(final_output, logger=logger)
        working_output.replace(final_output)

    cleanup_diagnostics(
        sidecars.values(),
        diagnostics_enabled=diagnostics_enabled,
        working_root=working_dir,
        final_root=final_dir,
        canonical_prefix=canonical_prefix,
        canonical_allowed=canonical_allowed,
        logger=logger,
    )

    for sidecar in sidecars.values():
        working_path = sidecar.working_path
        destination = sidecar.destination
        if working_path is None or not working_path.exists():
            continue
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
                remove_path(final_path, logger=logger)
                cleanup_empty_directories(final_parent, root=final_dir)
            continue
        if destination.exists():
            remove_path(destination, logger=logger)
        original_parent = working_path.parent
        working_path.replace(destination)
        cleanup_empty_directories(original_parent, root=working_dir)
        final_path = sidecar.final_path
        if final_path is not None and final_path.exists() and final_path != destination:
            final_parent = final_path.parent
            remove_path(final_path, logger=logger)
            cleanup_empty_directories(final_parent, root=final_dir)
    if sentinel_path.exists():
        remove_path(sentinel_path, logger=logger)


def resolve_postprocess_table(step: Any, final_output: Path) -> str | None:
    table = getattr(step, "output_stem", "")
    if not table or table not in POSTPROCESS_SUPPORTED_TABLES:
        return None
    if final_output.suffix.lower() != ".csv":
        return None
    return str(table)


def run_postprocess_hook(
    step: Any,
    final_output: Path,
    *,
    table: str | None = None,
    logger: logging.Logger,
) -> PostprocessResult | None:
    """Execute the post-processing pipeline for ``step`` when available."""

    if table is None:
        table = resolve_postprocess_table(step, final_output)
    if table is None:
        return None
    if not final_output.exists():
        logger.warning(
            "postprocess_input_missing",
            step=getattr(step, "name", "unknown"),
            table=table,
            input=str(final_output),
        )
        return None

    handlers = _POSTPROCESS_HANDLERS.get(table)
    if handlers is None:
        logger.warning(
            "postprocess_handlers_missing",
            step=getattr(step, "name", "unknown"),
            table=table,
        )
        return None
    destination = final_output.with_name(f"output_postprocessed.{table}.csv")
    logger.info(
        "postprocess_start",
        step=getattr(step, "name", "unknown"),
        table=table,
        input=str(final_output),
        output=str(destination),
    )
    pipeline_cfg = load_postprocess_pipeline_config(table, None)
    csv_cfg = get_postprocess_csv_config(pipeline_cfg)
    runtime_config = PostprocessingPipelineConfig(
        pipeline_config=pipeline_cfg,
        csv_runtime_config=csv_cfg,
        runner=handlers.runner,
        validator=handlers.validator,
        schema=handlers.schema,
        logger=logger,
    )
    pipeline_result = run_postprocessing_pipeline(
        table,
        final_output,
        destination,
        runtime_config,
    )
    metrics = pipeline_result.metrics
    report_path = pipeline_result.report_path
    if metrics is None:
        raise RuntimeError(f"postprocess metrics missing for table {table}")
    if report_path is None:
        raise RuntimeError(f"postprocess report missing for table {table}")
    logger.info(
        "postprocess_done",
        step=getattr(step, "name", "unknown"),
        table=table,
        output=str(pipeline_result.output_path),
        report=str(report_path),
        rows=int(metrics.output_rows),
        columns=int(metrics.output_columns),
    )
    return PostprocessResult(
        table=table,
        output_path=pipeline_result.output_path,
        report_path=report_path,
        metrics=metrics,
    )


@dataclass(frozen=True)
class _PostprocessHandlers:
    runner: Any
    validator: Any
    schema: Any


_POSTPROCESS_HANDLERS: dict[str, _PostprocessHandlers] = {
    "documents": _PostprocessHandlers(
        runner=run_document_pipeline,
        validator=validate_document_reports,
        schema=DOCUMENT_SCHEMA,
    ),
    "targets": _PostprocessHandlers(
        runner=run_target_pipeline,
        validator=validate_target_reports,
        schema=TARGET_SCHEMA,
    ),
    "assays": _PostprocessHandlers(
        runner=run_assay_pipeline,
        validator=validate_assay_reports,
        schema=ASSAY_SCHEMA,
    ),
    "activities": _PostprocessHandlers(
        runner=run_activity_pipeline,
        validator=validate_activity_reports,
        schema=ACTIVITY_SCHEMA,
    ),
    "testitem": _PostprocessHandlers(
        runner=run_testitem_pipeline,
        validator=validate_testitems,
        schema=TESTITEM_SCHEMA,
    ),
    "testitems": _PostprocessHandlers(
        runner=run_testitem_pipeline,
        validator=validate_testitems,
        schema=TESTITEM_SCHEMA,
    ),
}

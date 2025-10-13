"""Utilities for pruning legacy artefacts after the retention policy change."""

from __future__ import annotations

import json
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING

from ..common.log import logger as default_logger
from ..project_version import get_pipeline_version

if TYPE_CHECKING:  # pragma: no cover - imported for typing only
    from ..cli import Logger
    from ..config import Config


DEFAULT_PATTERNS: tuple[str, ...] = (
    "*.meta.yaml",
    "*.meta.yaml.lock",
    "*.quality.json",
    "*.quality.json.lock",
    "*_failure_cases.csv",
    "*_failure_cases.csv.meta.yaml",
    "*_failure_cases.csv.lock",
)

SENTINEL_FILENAME = ".chembl-da-legacy-cleanup.json"


@dataclass(frozen=True, slots=True)
class CleanupResult:
    """Summary of a legacy artefact cleanup run."""

    output_dir: Path
    sentinel_path: Path
    removed: tuple[Path, ...]
    errors: tuple[tuple[Path, str], ...]
    skipped: bool
    dry_run: bool

    @property
    def performed(self) -> bool:
        """Return ``True`` when cleanup attempted (even if nothing removed)."""

        return not self.skipped

    @property
    def removed_count(self) -> int:
        """Number of artefacts deleted (or slated for deletion)."""

        return len(self.removed)

    @property
    def error_count(self) -> int:
        """Number of artefacts that could not be removed."""

        return len(self.errors)


def list_legacy_artifacts(
    output_dir: Path, *, patterns: Sequence[str] | None = None
) -> list[Path]:
    """Return sorted artefacts within ``output_dir`` matching ``patterns``."""

    resolved_dir = Path(output_dir)
    if not resolved_dir.exists():
        return []

    use_patterns = tuple(patterns) if patterns else DEFAULT_PATTERNS
    seen: set[Path] = set()
    for pattern in use_patterns:
        for candidate in resolved_dir.rglob(pattern):
            if candidate.is_file():
                seen.add(candidate.resolve())
    return sorted(seen)


def _write_sentinel(
    path: Path, *, removed: Iterable[Path], errors: Iterable[tuple[Path, str]]
) -> None:
    data = {
        "version": get_pipeline_version(),
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "removed": [str(item) for item in removed],
        "errors": [{"path": str(path), "error": message} for path, message in errors],
    }
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def cleanup_legacy_outputs(
    output_dir: Path,
    *,
    patterns: Sequence[str] | None = None,
    dry_run: bool = False,
    force: bool = False,
    logger: Logger | None = None,
) -> CleanupResult:
    """Remove legacy artefacts from ``output_dir`` once per installation."""

    use_logger = logger or default_logger
    resolved_dir = Path(output_dir)
    sentinel = resolved_dir / SENTINEL_FILENAME

    if not resolved_dir.exists():
        use_logger.debug(
            "legacy_cleanup_skipped",
            reason="output_dir_missing",
            directory=str(resolved_dir),
        )
        return CleanupResult(
            output_dir=resolved_dir,
            sentinel_path=sentinel,
            removed=(),
            errors=(),
            skipped=True,
            dry_run=dry_run,
        )

    if sentinel.exists() and not force:
        use_logger.debug(
            "legacy_cleanup_skipped",
            reason="sentinel_present",
            directory=str(resolved_dir),
            sentinel=str(sentinel),
        )
        return CleanupResult(
            output_dir=resolved_dir,
            sentinel_path=sentinel,
            removed=(),
            errors=(),
            skipped=True,
            dry_run=dry_run,
        )

    artefacts = list_legacy_artifacts(resolved_dir, patterns=patterns)
    removed: list[Path] = []
    errors: list[tuple[Path, str]] = []

    for candidate in artefacts:
        if dry_run:
            removed.append(candidate)
            continue
        try:
            Path(candidate).unlink(missing_ok=False)
        except FileNotFoundError:
            continue
        except OSError as exc:  # pragma: no cover - rare, environment specific
            errors.append((Path(candidate), str(exc)))
            use_logger.warning(
                "legacy_cleanup_error",
                path=str(candidate),
                error=str(exc),
            )
        else:
            removed.append(Path(candidate))

    if not dry_run:
        try:
            _write_sentinel(sentinel, removed=removed, errors=errors)
        except OSError as exc:  # pragma: no cover - disk error
            use_logger.warning(
                "legacy_cleanup_sentinel_failed",
                path=str(sentinel),
                error=str(exc),
            )

    if removed:
        use_logger.info(
            "legacy_cleanup_removed",
            removed=len(removed),
            directory=str(resolved_dir),
            dry_run=dry_run,
        )
    else:
        use_logger.debug(
            "legacy_cleanup_noop",
            directory=str(resolved_dir),
            dry_run=dry_run,
        )

    return CleanupResult(
        output_dir=resolved_dir,
        sentinel_path=sentinel,
        removed=tuple(removed),
        errors=tuple(errors),
        skipped=False,
        dry_run=dry_run,
    )


def ensure_legacy_cleanup(
    cfg: Config, *, logger: Logger | None = None, dry_run: bool = False
) -> CleanupResult:
    """Ensure the legacy cleanup runs once using the configured output dir."""

    output_dir = Path(cfg.io.output_dir)
    return cleanup_legacy_outputs(
        output_dir,
        dry_run=dry_run,
        logger=logger,
    )


__all__ = [
    "CleanupResult",
    "DEFAULT_PATTERNS",
    "SENTINEL_FILENAME",
    "cleanup_legacy_outputs",
    "ensure_legacy_cleanup",
    "list_legacy_artifacts",
]

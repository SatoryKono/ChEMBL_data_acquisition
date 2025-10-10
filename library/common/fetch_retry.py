"""Helpers for retrying chunked fetch operations and tracking failures."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ..config import Config, RetryCfg
from ..sidecar import SidecarErrors

_CHUNK_FAILURE_IDS_LIMIT = 100


@dataclass(slots=True)
class _ChunkFailure:
    """Internal representation of a failed chunk fetch."""

    chunk_ids: list[str]
    error: str


class ChunkFailureTracker:
    """Collect chunk level failures and expose aggregated statistics."""

    def __init__(self) -> None:
        self._failures: list[_ChunkFailure] = []
        self._sidecar = SidecarErrors()

    def add_failure(self, chunk_ids: Sequence[str], error: str) -> None:
        """Record a failed chunk fetch attempt."""

        ids = [str(identifier) for identifier in chunk_ids]
        failure = _ChunkFailure(chunk_ids=ids, error=error)
        self._failures.append(failure)
        self._sidecar.add_error(
            {
                "chunk_ids": ",".join(ids),
                "chunk_size": len(ids),
                "error": error,
            }
        )

    def stats(self) -> dict[str, Any]:
        """Return aggregated statistics for persisted metadata files."""

        if not self._failures:
            return {}
        unique_ids = sorted(
            {
                identifier
                for failure in self._failures
                for identifier in failure.chunk_ids
            }
        )
        total_unique_ids = len(unique_ids)
        truncated = total_unique_ids > _CHUNK_FAILURE_IDS_LIMIT
        if truncated:
            reported_ids = list(unique_ids[:_CHUNK_FAILURE_IDS_LIMIT])
        else:
            reported_ids = list(unique_ids)
        return {
            "chunk_fetch_failure_chunks": len(self._failures),
            "chunk_fetch_failure_ids": reported_ids,
            "chunk_fetch_failure_ids_total": total_unique_ids,
            "chunk_fetch_failure_ids_truncated": truncated,
        }

    def save(self, path: Path, *, cfg: Config | None = None) -> None:
        """Write recorded failures to ``path`` and emit metadata sidecar."""

        if self._failures:
            self._sidecar.save(path, cfg=cfg)
            return
        path.unlink(missing_ok=True)
        Path(f"{path}.meta.yaml").unlink(missing_ok=True)

    def has_failures(self) -> bool:
        """Return ``True`` when at least one failure has been recorded."""

        return bool(self._failures)


def compute_backoff_delay(
    attempt: int,
    retry_cfg: RetryCfg,
    *,
    jitter: Callable[[float], float] | None = None,
) -> float:
    """Return the exponential backoff delay for ``attempt`` respecting caps."""

    if attempt <= 0:
        raise ValueError("attempt must be a positive integer")
    factor = retry_cfg.backoff_factor
    base_delay = factor * (2 ** (attempt - 1))
    delay = base_delay
    if jitter is not None:
        delay += jitter(base_delay)
    cap = retry_cfg.backoff_cap
    if cap is not None:
        delay = min(delay, cap)
    return delay


__all__ = ["ChunkFailureTracker", "compute_backoff_delay"]

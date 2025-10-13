"""Helpers for retrying chunked fetch operations and tracking failures."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator, Sequence
from dataclasses import dataclass
from itertools import islice
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
    sidecar_row: dict[str, Any]


class ChunkFailureTracker:
    """Collect chunk level failures and expose aggregated statistics."""

    def __init__(self) -> None:
        self._failures: list[_ChunkFailure] = []
        self._sidecar = SidecarErrors()

    def add_failure(self, chunk_ids: Sequence[str], error: str) -> None:
        """Record a failed chunk fetch attempt."""

        ids = [str(identifier) for identifier in chunk_ids]
        sidecar_row = {
            "chunk_ids": ",".join(ids),
            "chunk_size": len(ids),
            "error": error,
            "chunk_fetch_failure_ids_total": 0,
            "chunk_fetch_failure_ids_truncated": False,
        }
        self._sidecar.add_error(sidecar_row)
        failure = _ChunkFailure(chunk_ids=ids, error=error, sidecar_row=sidecar_row)
        self._failures.append(failure)

    def stats(self) -> dict[str, Any]:
        """Return aggregated statistics for persisted metadata files.

        Only the first ``_CHUNK_FAILURE_IDS_LIMIT`` unique identifiers are
        materialised to keep the memory footprint bounded. The full unique count
        is still computed to preserve the historical statistics exported in
        sidecars and logs.
        """

        if not self._failures:
            return {}

        unique_ids_iter = _iter_unique_ids(_iter_all_chunk_ids(self._failures))
        reported_ids = list(islice(unique_ids_iter, _CHUNK_FAILURE_IDS_LIMIT))
        total_unique_ids = len(reported_ids) + sum(1 for _ in unique_ids_iter)
        truncated = total_unique_ids > len(reported_ids)

        self._update_sidecar_summary(total_unique_ids, truncated)
        return {
            "chunk_fetch_failure_chunks": len(self._failures),
            "chunk_fetch_failure_ids": reported_ids,
            "chunk_fetch_failure_ids_total": total_unique_ids,
            "chunk_fetch_failure_ids_truncated": truncated,
        }

    def save(self, path: Path, *, cfg: Config | None = None) -> None:
        """Write recorded failures to ``path`` and emit metadata sidecar."""

        if self._failures:
            self.stats()
            self._sidecar.save(path, cfg=cfg)
            return
        path.unlink(missing_ok=True)
        Path(f"{path}.meta.yaml").unlink(missing_ok=True)

    def has_failures(self) -> bool:
        """Return ``True`` when at least one failure has been recorded."""

        return bool(self._failures)

    def _update_sidecar_summary(self, total: int, truncated: bool) -> None:
        for failure in self._failures:
            failure.sidecar_row["chunk_fetch_failure_ids_total"] = total
            failure.sidecar_row["chunk_fetch_failure_ids_truncated"] = truncated


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


def _iter_all_chunk_ids(failures: Sequence[_ChunkFailure]) -> Iterator[str]:
    for failure in failures:
        yield from failure.chunk_ids


def _iter_unique_ids(identifiers: Iterable[str]) -> Iterator[str]:
    seen: set[str] = set()
    for identifier in identifiers:
        if identifier in seen:
            continue
        seen.add(identifier)
        yield identifier

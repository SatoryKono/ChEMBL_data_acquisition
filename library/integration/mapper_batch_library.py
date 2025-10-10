"""Batch mapping utilities for ChEMBL to UniProt conversion."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import TypeVar

from ..common.log import logger
from ..common.rate_limiter import get_limiter
from ..config import UniprotMappingCfg
from .mapper_library import map_chembl_to_uniprot

T = TypeVar("T")


def batch_iterable(iterable: Iterable[T], batch_size: int) -> Iterator[list[T]]:
    """Yield items from ``iterable`` in batches of ``batch_size``.

    Parameters
    ----------
    iterable:
        Source of items to batch.
    batch_size:
        Maximum size of each batch. Must be greater than zero.

    Yields
    ------
    list
        Lists containing up to ``batch_size`` elements from ``iterable``.
    """
    if batch_size <= 0:
        raise ValueError("batch_size must be greater than zero")

    batch: list[T] = []
    for item in iterable:
        batch.append(item)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def map_chembl_ids_to_uniprot(
    chembl_ids: Iterable[str],
    cfg: UniprotMappingCfg,
    *,
    batch_size: int,
    rps: float,
    max_workers: int | None = None,
) -> dict[str, str | None]:
    """Map multiple ChEMBL identifiers to UniProt accessions.

    Parameters
    ----------
    chembl_ids:
        Iterable of ChEMBL target identifiers.
    cfg:
        Configuration for the UniProt mapping API.
    batch_size:
        Number of identifiers processed per worker.
    rps:
        Maximum requests per second across all workers.
    max_workers:
        Optional upper bound on the number of worker threads. Defaults to
        ``None`` which lets :class:`~concurrent.futures.ThreadPoolExecutor`
        decide.

    Returns
    -------
    dict
        Mapping of input ChEMBL IDs to UniProt accessions. If mapping fails for
        an ID the value will be ``None``.
    """
    limiter = get_limiter("chembl_uniprot_batch", rps)

    def _process(batch: list[str]) -> dict[str, str | None]:
        results: dict[str, str | None] = {}
        for chembl_id in batch:
            limiter.acquire()
            try:
                results[chembl_id] = map_chembl_to_uniprot(chembl_id, cfg)
            except Exception as exc:  # pragma: no cover - network issues
                logger.warning("batch_map_failed", chembl_id=chembl_id, error=str(exc))
                results[chembl_id] = None
        return results

    results: dict[str, str | None] = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(_process, b) for b in batch_iterable(chembl_ids, batch_size)
        ]
        for future in as_completed(futures):
            batch_result = future.result()
            results.update(batch_result)
    return results

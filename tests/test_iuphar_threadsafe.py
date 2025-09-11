"""Thread-safety tests for :mod:`library.iuphar_library`."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import threading

import pandas as pd

import library.iuphar_library as ii
from library.config import IupharCfg


class DummyResponse:
    def __enter__(self) -> "DummyResponse":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # pragma: no cover - no errors
        return None

    def raise_for_status(self) -> None:
        return None

    def json(self) -> list[dict[str, int]]:
        return [{"id": 1}]


class DummySession:
    def __init__(self) -> None:
        self.calls: list[int] = []

    def get(self, url: str, timeout: tuple[int, int]) -> DummyResponse:
        self.calls.append(threading.get_ident())
        return DummyResponse()


class DummyLimiter:
    def acquire(self) -> None:  # pragma: no cover - no delay
        return None


def test_websearch_gene_to_id_thread_safe(monkeypatch) -> None:
    data = ii.IUPHARData(target_df=pd.DataFrame(), family_df=pd.DataFrame())
    dummy = DummySession()
    monkeypatch.setattr(ii, "_session", dummy)
    monkeypatch.setattr(ii, "get_limiter", lambda name, rps: DummyLimiter())

    cfg = IupharCfg(base="https://example.org/services", rps=10, burst=10)

    def worker() -> dict:
        return data.websearch_gene_to_id("GENE", cfg)

    with ThreadPoolExecutor(max_workers=10) as pool:
        results = [f.result() for f in (pool.submit(worker) for _ in range(10))]

    assert all(r == {"id": 1} for r in results)
    assert len(dummy.calls) == 10

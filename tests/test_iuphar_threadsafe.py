"""Thread-safety tests for :mod:`library.iuphar_library`."""

from __future__ import annotations

import threading
import time
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import pytest

from library import iuphar_library as ii
from library.clients import iuphar as ci
from library.config import ApiCfg, IupharCfg, RetryCfg


class DummyResponse:
    def __enter__(self) -> DummyResponse:
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
    monkeypatch.setattr(ci, "_session", dummy)
    monkeypatch.setattr(ci, "get_limiter", lambda name, rps, burst=None: DummyLimiter())

    cfg = IupharCfg(base="https://example.org/services", rps=10, burst=10)

    def worker() -> dict:
        return data.websearch_gene_to_id("GENE", cfg)

    with ThreadPoolExecutor(max_workers=10) as pool:
        results = [f.result() for f in (pool.submit(worker) for _ in range(10))]

    assert all(r == {"id": 1} for r in results)
    assert len(dummy.calls) == 10


"""Thread-safety tests for IUPHAR session handling."""


def test_session_serialization(monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure concurrent calls to the shared session are serialised."""
    data = ii.IUPHARData(target_df=pd.DataFrame(), family_df=pd.DataFrame())
    cfg = IupharCfg(
        base="https://example.org", timeout_connect=1, timeout_read=1, rps=1, burst=1
    )

    class DummyLimiter:
        def acquire(self) -> None:  # pragma: no cover - trivial
            return None

    monkeypatch.setattr(ci, "get_limiter", lambda *a, **k: DummyLimiter())

    order: list[str] = []
    in_use = False
    lock = threading.Lock()

    def fake_get(url: str, timeout: tuple[int, int]):
        nonlocal in_use
        with lock:
            assert not in_use, "Concurrent session usage detected"
            in_use = True
            order.append("start")
        time.sleep(0.05)
        with lock:
            in_use = False
            order.append("end")

        class Resp:
            def __enter__(self):
                return self

            def __exit__(self, *exc):
                return False

            def raise_for_status(self) -> None:
                return None

            def json(self):  # pragma: no cover - unused
                return []

        return Resp()

    monkeypatch.setattr(ci._session, "get", fake_get)

    threads = [
        threading.Thread(target=data.websearch_gene_to_id, args=("GENE", cfg))
        for _ in range(2)
    ]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert order == ["start", "end", "start", "end"]


def test_init_session_closes_previous_session(monkeypatch: pytest.MonkeyPatch) -> None:
    """Reinitialising the session should close the previous instance exactly once."""

    closed: list[str] = []

    class DummySession:
        def __init__(self, name: str) -> None:
            self.name = name
            self.closed = False

        def close(self) -> None:
            self.closed = True
            closed.append(self.name)

    dummy_initial = DummySession("initial")
    monkeypatch.setattr(ci, "_session", dummy_initial, raising=False)

    def fake_session_with_retry(api: ApiCfg, retry: RetryCfg) -> DummySession:
        return DummySession(api.user_agent)

    monkeypatch.setattr(ci, "session_with_retry", fake_session_with_retry)

    retry_cfg = RetryCfg()
    first_api = ApiCfg(user_agent="chembl-tests/1.0 (mailto:tests@ebi.ac.uk)")
    second_api = ApiCfg(user_agent="chembl-tests/2.0 (mailto:tests@ebi.ac.uk)")

    ci.init_session(first_api, retry_cfg)
    first_session = ci._session

    assert isinstance(first_session, DummySession)
    assert closed == ["initial"]
    assert dummy_initial.closed is True

    ci.init_session(second_api, retry_cfg)
    second_session = ci._session

    assert isinstance(second_session, DummySession)
    assert closed == ["initial", first_api.user_agent]
    assert first_session.closed is True
    assert second_session.closed is False

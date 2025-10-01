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
        self.owner: int | None = None

    def get(self, url: str, timeout: tuple[int, int]) -> DummyResponse:
        thread_id = threading.get_ident()
        if self.owner is None:
            self.owner = thread_id
        else:
            assert (
                self.owner == thread_id
            ), "Session reused across threads"
        self.calls.append(thread_id)
        return DummyResponse()


class DummyLimiter:
    def acquire(self) -> None:  # pragma: no cover - no delay
        return None


def test_websearch_gene_to_id_thread_safe(monkeypatch) -> None:
    data = ii.IUPHARData(target_df=pd.DataFrame(), family_df=pd.DataFrame())
    call_threads: list[int] = []
    created_sessions: list[DummySession] = []

    def session_factory() -> DummySession:
        session = DummySession()
        created_sessions.append(session)
        return session

    monkeypatch.setattr(ci, "_session_factory", session_factory, raising=False)
    monkeypatch.setattr(ci, "_session_local", threading.local(), raising=False)
    monkeypatch.setattr(ci, "_sessions", set(), raising=False)
    monkeypatch.setattr(ci, "_active_requests", 0, raising=False)
    monkeypatch.setattr(
        ci, "get_limiter", lambda name, rps, burst=None: DummyLimiter()
    )

    cfg = IupharCfg(base="https://example.org/services", rps=10, burst=10)
    genes = [f"GENE-{idx}" for idx in range(10)]

    def worker(symbol: str) -> dict:
        result = data.websearch_gene_to_id(symbol, cfg)
        call_threads.append(threading.get_ident())
        return result

    with ThreadPoolExecutor(max_workers=10) as pool:
        futures = [pool.submit(worker, gene) for gene in genes]
        results = [f.result() for f in futures]

    assert all(r == {"id": 1} for r in results)
    assert len(call_threads) == 10
    assert created_sessions
    owners = [session.owner for session in created_sessions if session.owner is not None]
    assert len(set(owners)) == len(owners)
    total_calls = sum(len(session.calls) for session in created_sessions)
    assert total_calls == len(call_threads)


"""Thread-safety tests for IUPHAR session handling."""


def test_session_is_thread_local(monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure each thread obtains an independent session instance."""

    data = ii.IUPHARData(target_df=pd.DataFrame(), family_df=pd.DataFrame())
    cfg = IupharCfg(
        base="https://example.org", timeout_connect=1, timeout_read=1, rps=1, burst=1
    )

    class DummyLimiter:
        def acquire(self) -> None:  # pragma: no cover - trivial
            return None

    monkeypatch.setattr(ci, "get_limiter", lambda *a, **k: DummyLimiter())

    created_sessions: list[DummySession] = []

    def session_factory() -> DummySession:
        session = DummySession()
        created_sessions.append(session)
        return session

    monkeypatch.setattr(ci, "_session_factory", session_factory, raising=False)
    monkeypatch.setattr(ci, "_session_local", threading.local(), raising=False)
    monkeypatch.setattr(ci, "_sessions", set(), raising=False)
    monkeypatch.setattr(ci, "_active_requests", 0, raising=False)

    barrier = threading.Barrier(2)

    def worker(symbol: str) -> dict:
        barrier.wait()
        return data.websearch_gene_to_id(symbol, cfg)

    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(worker, symbol) for symbol in ("GENE-A", "GENE-B")]
        [f.result() for f in futures]

    owners = {session.owner for session in created_sessions if session.owner is not None}

    assert len(created_sessions) == 2
    assert len(owners) == 2


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
    monkeypatch.setattr(ci, "_session_factory", lambda: dummy_initial, raising=False)
    monkeypatch.setattr(ci, "_session_local", threading.local(), raising=False)
    monkeypatch.setattr(ci, "_sessions", {dummy_initial}, raising=False)
    monkeypatch.setattr(ci, "_active_requests", 0, raising=False)

    def fake_session_with_retry(api: ApiCfg, retry: RetryCfg) -> DummySession:
        return DummySession(api.user_agent)

    monkeypatch.setattr(ci, "session_with_retry", fake_session_with_retry)

    retry_cfg = RetryCfg()
    first_api = ApiCfg(user_agent="chembl-tests/1.0 (mailto:tests@ebi.ac.uk)")
    second_api = ApiCfg(user_agent="chembl-tests/2.0 (mailto:tests@ebi.ac.uk)")

    ci.init_session(first_api, retry_cfg)
    first_session = ci.get_session()

    assert isinstance(first_session, DummySession)
    assert closed == ["initial"]
    assert dummy_initial.closed is True

    ci.init_session(second_api, retry_cfg)
    second_session = ci.get_session()

    assert isinstance(second_session, DummySession)
    assert closed == ["initial", first_api.user_agent]
    assert first_session.closed is True
    assert second_session.closed is False

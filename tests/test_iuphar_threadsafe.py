"""Thread-safety tests for IUPHAR session handling."""

from __future__ import annotations

import threading
import time

import pandas as pd
import pytest

from library import iuphar_library as ii
from library.config import IupharCfg


def test_session_serialization(monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure concurrent calls to the shared session are serialised."""
    data = ii.IUPHARData(target_df=pd.DataFrame(), family_df=pd.DataFrame())
    cfg = IupharCfg(
        base="https://example.org", timeout_connect=1, timeout_read=1, rps=0, burst=1
    )

    class DummyLimiter:
        def acquire(self) -> None:  # pragma: no cover - trivial
            return None

    monkeypatch.setattr(ii, "get_limiter", lambda *a, **k: DummyLimiter())

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

    monkeypatch.setattr(ii._session, "get", fake_get)

    threads = [
        threading.Thread(target=data.websearch_gene_to_id, args=("GENE", cfg))
        for _ in range(2)
    ]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert order == ["start", "end", "start", "end"]

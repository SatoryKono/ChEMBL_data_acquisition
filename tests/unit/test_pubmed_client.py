from __future__ import annotations

import pytest
from urllib3.util import retry as urllib3_retry

from library.clients import pubmed
from library.config import ApiCfg, RetryCfg, session_with_retry


@pytest.mark.unit
def test_retry_delay__respects_backoff_cap(monkeypatch: pytest.MonkeyPatch) -> None:
    retry_cfg = RetryCfg(max_attempts=4, backoff_factor=2.0, backoff_cap=3.0)

    monkeypatch.setattr(pubmed.random, "uniform", lambda *args, **kwargs: 0.5)

    delay = pubmed._retry_delay(2, 0.25, retry_cfg, timeout=None)

    assert delay == pytest.approx(3.0)


@pytest.mark.unit
def test_session_with_retry__uses_backoff_cap() -> None:
    api_cfg = ApiCfg()
    retry_cfg = RetryCfg(max_attempts=4, backoff_factor=1.0, backoff_cap=7.5)

    with session_with_retry(api_cfg, retry_cfg) as session:
        adapter = session.get_adapter("https://")
        max_retries = adapter.max_retries

    assert isinstance(max_retries, urllib3_retry.Retry)
    assert max_retries.backoff_max == pytest.approx(7.5)

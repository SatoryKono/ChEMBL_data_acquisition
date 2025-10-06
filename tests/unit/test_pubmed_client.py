from __future__ import annotations

import pytest
from urllib3.util import retry as urllib3_retry

from library.clients import pubmed
from library.config import ApiCfg, RetryCfg, session_with_retry


@pytest.mark.unit
def test_retry_delay__respects_backoff_cap() -> None:
    retry_cfg = RetryCfg(max_attempts=4, backoff_factor=2.0, backoff_cap=3.0)

    delay = pubmed._retry_delay(
        2,
        0.25,
        retry_cfg,
        timeout=None,
        jitter=lambda _: 0.5,
    )

    assert delay == pytest.approx(3.0)


@pytest.mark.unit
def test_retry_delay__deterministic_jitter() -> None:
    retry_cfg = RetryCfg(max_attempts=4, backoff_factor=1.0, backoff_cap=None, jitter_seed=7)
    jitter_one = retry_cfg.build_jitter()
    jitter_two = RetryCfg(max_attempts=4, backoff_factor=1.0, backoff_cap=None, jitter_seed=7).build_jitter()

    assert jitter_one is not None
    assert jitter_two is not None

    base_delay = 0.25
    delays_one = [
        pubmed._retry_delay(attempt, base_delay, retry_cfg, timeout=None, jitter=jitter_one)
        for attempt in range(1, 4)
    ]
    delays_two = [
        pubmed._retry_delay(attempt, base_delay, retry_cfg, timeout=None, jitter=jitter_two)
        for attempt in range(1, 4)
    ]

    assert delays_one == delays_two


@pytest.mark.unit
def test_session_with_retry__uses_backoff_cap() -> None:
    api_cfg = ApiCfg()
    retry_cfg = RetryCfg(max_attempts=4, backoff_factor=1.0, backoff_cap=7.5)

    with session_with_retry(api_cfg, retry_cfg) as session:
        adapter = session.get_adapter("https://")
        max_retries = adapter.max_retries

    assert isinstance(max_retries, urllib3_retry.Retry)
    assert max_retries.backoff_max == pytest.approx(7.5)

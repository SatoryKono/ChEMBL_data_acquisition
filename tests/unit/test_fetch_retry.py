from __future__ import annotations

import pytest

from library.common.fetch_retry import compute_backoff_delay
from library.config import RetryCfg


def test_compute_backoff_delay__deterministic_with_seed() -> None:
    retry_cfg = RetryCfg(backoff_factor=1.25, backoff_cap=None, jitter_seed=7)
    delays_first = [compute_backoff_delay(attempt, retry_cfg) for attempt in range(1, 5)]

    retry_cfg_same_seed = RetryCfg(backoff_factor=1.25, backoff_cap=None, jitter_seed=7)
    delays_second = [compute_backoff_delay(attempt, retry_cfg_same_seed) for attempt in range(1, 5)]

    assert delays_first == delays_second


def test_compute_backoff_delay__adds_jitter_before_cap() -> None:
    retry_cfg = RetryCfg(backoff_factor=2.0, backoff_cap=3.5, jitter_seed=11)
    jitter = RetryCfg(backoff_factor=2.0, backoff_cap=3.5, jitter_seed=11).build_jitter()
    assert jitter is not None

    attempt = 2
    base_delay = retry_cfg.backoff_factor * (2 ** (attempt - 1))
    expected = min(base_delay + jitter(retry_cfg.backoff_factor), retry_cfg.backoff_cap)

    assert compute_backoff_delay(attempt, retry_cfg) == pytest.approx(expected)

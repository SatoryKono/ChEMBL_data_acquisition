from __future__ import annotations

import pytest

from library.common.fetch_retry import ChunkFailureTracker, compute_backoff_delay
from library.config import RetryCfg


def test_compute_backoff_delay__deterministic_with_cached_jitter() -> None:
    retry_cfg = RetryCfg(backoff_factor=1.25, backoff_cap=None, jitter_seed=7)
    jitter = retry_cfg.build_jitter()
    assert jitter is not None
    delays_first = [
        compute_backoff_delay(attempt, retry_cfg, jitter=jitter)
        for attempt in range(1, 5)
    ]

    retry_cfg_same_seed = RetryCfg(backoff_factor=1.25, backoff_cap=None, jitter_seed=7)
    jitter_same_seed = retry_cfg_same_seed.build_jitter()
    assert jitter_same_seed is not None
    delays_second = [
        compute_backoff_delay(attempt, retry_cfg_same_seed, jitter=jitter_same_seed)
        for attempt in range(1, 5)
    ]

    assert delays_first == delays_second


def test_compute_backoff_delay__jitter_sequence_restarts_with_same_seed() -> None:
    retry_cfg = RetryCfg(backoff_factor=1.0, backoff_cap=None, jitter_seed=17)
    jitter = retry_cfg.build_jitter()
    assert jitter is not None

    attempts = range(1, 5)
    first_run = [
        compute_backoff_delay(attempt, retry_cfg, jitter=jitter) for attempt in attempts
    ]

    base_delays = [
        retry_cfg.backoff_factor * (2 ** (attempt - 1)) for attempt in attempts
    ]
    jitter_offsets = [
        value - base for value, base in zip(first_run, base_delays, strict=True)
    ]

    assert jitter_offsets[0] != pytest.approx(jitter_offsets[1])

    jitter_repeat = retry_cfg.build_jitter()
    assert jitter_repeat is not None
    second_run = [
        compute_backoff_delay(attempt, retry_cfg, jitter=jitter_repeat)
        for attempt in attempts
    ]

    assert first_run == second_run


def test_compute_backoff_delay__adds_jitter_before_cap() -> None:
    retry_cfg = RetryCfg(backoff_factor=2.0, backoff_cap=3.5, jitter_seed=11)
    jitter = retry_cfg.build_jitter()
    assert jitter is not None
    jitter_expected = retry_cfg.build_jitter()
    assert jitter_expected is not None

    attempt = 2
    base_delay = retry_cfg.backoff_factor * (2 ** (attempt - 1))
    expected = min(
        base_delay + jitter_expected(retry_cfg.backoff_factor),
        retry_cfg.backoff_cap,
    )

    assert compute_backoff_delay(attempt, retry_cfg, jitter=jitter) == pytest.approx(
        expected
    )


def test_chunk_failure_tracker_stats__returns_fresh_mapping_when_empty() -> None:
    tracker = ChunkFailureTracker()

    first = tracker.stats()
    first["custom"] = "value"

    second = tracker.stats()

    assert first is not second
    assert second == {}


def test_chunk_failure_tracker_stats__does_not_share_lists() -> None:
    tracker = ChunkFailureTracker()
    tracker.add_failure(["A", "B"], "boom")

    first = tracker.stats()
    first["chunk_fetch_failure_ids"].append("C")

    second = tracker.stats()

    assert second["chunk_fetch_failure_ids"] == ["A", "B"]
    assert second["chunk_fetch_failure_ids_total"] == 2
    assert second["chunk_fetch_failure_ids_truncated"] is False


def test_chunk_failure_tracker_stats__limits_reported_ids() -> None:
    tracker = ChunkFailureTracker()

    for index in range(150):
        tracker.add_failure([f"ID{index:03d}"], "boom")

    stats = tracker.stats()

    assert stats["chunk_fetch_failure_ids_total"] == 150
    assert stats["chunk_fetch_failure_ids_truncated"] is True
    assert (
        len(stats["chunk_fetch_failure_ids"]) < stats["chunk_fetch_failure_ids_total"]
    )

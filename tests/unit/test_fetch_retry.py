from __future__ import annotations

import csv
import tracemalloc
from itertools import cycle, islice

import pytest
from hypothesis import HealthCheck, given, seed, settings, strategies as st

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
        base_delay + jitter_expected(base_delay),
        retry_cfg.backoff_cap,
    )

    assert compute_backoff_delay(attempt, retry_cfg, jitter=jitter) == pytest.approx(
        expected
    )


@pytest.mark.parametrize(
    ("attempt", "expected"),
    [
        (1, 1.0 + 0.25),
        (2, 2.0 + 0.5),
        (3, 4.0 + 1.0),
        (4, 8.0 + 2.0),
    ],
)
def test_compute_backoff_delay__growth_with_fixed_jitter(
    attempt: int, expected: float
) -> None:
    retry_cfg = RetryCfg(backoff_factor=1.0, backoff_cap=None)

    def _jitter(base_delay: float) -> float:
        return base_delay * 0.25

    assert compute_backoff_delay(attempt, retry_cfg, jitter=_jitter) == pytest.approx(
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


def test_chunk_failure_tracker_stats__preserves_first_hundred_ids_and_total() -> None:
    tracker = ChunkFailureTracker()

    all_identifiers = [f"ID{index:03d}" for index in range(150)]
    duplicated_stream: list[str] = []
    for position, identifier in enumerate(all_identifiers):
        duplicated_stream.append(identifier)
        if position % 10 == 0:
            duplicated_stream.append(identifier)

    chunk_size = 4
    for start in range(0, len(duplicated_stream), chunk_size):
        chunk = duplicated_stream[start : start + chunk_size]
        tracker.add_failure(chunk, "boom")

    stats = tracker.stats()

    assert stats["chunk_fetch_failure_ids"] == all_identifiers[:100]
    assert stats["chunk_fetch_failure_ids_total"] == len(all_identifiers)
    assert stats["chunk_fetch_failure_ids_truncated"] is True


@seed(42)
@settings(
    deadline=None,
    max_examples=3,
    suppress_health_check=[HealthCheck.too_slow],
)
@given(
    sample_size=st.integers(min_value=10_000, max_value=12_000),
    seed_values=st.lists(
        st.integers(min_value=0, max_value=50_000),
        min_size=256,
        max_size=512,
    ),
    chunk_size=st.integers(min_value=5, max_value=31),
)
def test_chunk_failure_tracker_stats__streaming_unique_ids_memory_bound(
    sample_size: int,
    seed_values: list[int],
    chunk_size: int,
) -> None:
    tracker = ChunkFailureTracker()

    prefix_length = sample_size // 2
    cycled = list(islice(cycle(seed_values), prefix_length))
    tail_start = 1_000_000
    tail = list(range(tail_start, tail_start + sample_size - prefix_length))
    formatted = [f"ID{value:05d}" for value in (cycled + tail)]
    for start in range(0, len(formatted), chunk_size):
        chunk = formatted[start : start + chunk_size]
        tracker.add_failure(chunk, "boom")

    tracemalloc.start()
    stats = tracker.stats()
    _current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    seen: set[str] = set()
    ordered_unique: list[str] = []
    for identifier in formatted:
        if identifier in seen:
            continue
        seen.add(identifier)
        ordered_unique.append(identifier)

    assert stats["chunk_fetch_failure_ids"] == ordered_unique[:100]
    assert stats["chunk_fetch_failure_ids_total"] == len(ordered_unique)
    assert stats["chunk_fetch_failure_ids_truncated"] is (len(ordered_unique) > 100)
    assert peak <= 5_000_000


def test_chunk_failure_tracker_save__sidecar_includes_truncation_metadata(
    tmp_path,
) -> None:
    tracker = ChunkFailureTracker()

    for index in range(120):
        tracker.add_failure([f"ID{index:03d}"], "boom")

    path = tmp_path / "failures.csv"
    tracker.save(path)

    with path.open("r", newline="", encoding="utf8") as fh:
        rows = list(csv.DictReader(fh))

    assert rows, "expected tracker.save to create sidecar entries"
    assert rows[0]["chunk_fetch_failure_ids_total"] == "120"
    assert rows[0]["chunk_fetch_failure_ids_truncated"] == "True"

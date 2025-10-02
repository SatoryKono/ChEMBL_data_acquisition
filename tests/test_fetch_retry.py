import csv
from pathlib import Path

import pytest

from library.config import RetryCfg
from library.fetch_retry import ChunkFailureTracker, compute_backoff_delay


def test_compute_backoff_delay_respects_cap() -> None:
    cfg = RetryCfg(max_attempts=3, backoff_factor=0.5, backoff_cap=1.0)
    assert compute_backoff_delay(1, cfg) == pytest.approx(0.5)
    assert compute_backoff_delay(3, cfg) == pytest.approx(1.0)


def test_compute_backoff_delay_raises_on_invalid_attempt() -> None:
    cfg = RetryCfg()
    with pytest.raises(ValueError):
        compute_backoff_delay(0, cfg)


def test_chunk_failure_tracker_records_stats(tmp_path: Path, cfg) -> None:
    tracker = ChunkFailureTracker()
    tracker.add_failure(["A1", "A2"], "boom")
    tracker.add_failure(["A2"], "retry failed")

    stats = tracker.stats()
    assert stats["chunk_fetch_failure_chunks"] == 2
    assert stats["chunk_fetch_failure_ids"] == ["A1", "A2"]

    path = tmp_path / "fetch_failures.csv"
    tracker.save(path, cfg=cfg)
    assert path.exists()
    with path.open(newline="", encoding="utf8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows and rows[0]["chunk_ids"] == "A1,A2"

    # Saving again without new failures should remove the artefact.
    empty_tracker = ChunkFailureTracker()
    path.write_text("stale", encoding="utf8")
    meta = Path(f"{path}.meta.yaml")
    meta.write_text("{}", encoding="utf8")
    empty_tracker.save(path, cfg=cfg)
    assert not path.exists()
    assert not meta.exists()

from __future__ import annotations

import io
import json
import sys

import pytest

from library.integration import mapper_batch_library as mbl
from library.cli import LoggerConfig, configure_logger
from library.config import UniprotMappingCfg


def test_batch_iterable() -> None:
    """Ensure items are grouped into correct batch sizes."""
    data = list(mbl.batch_iterable(range(5), 2))
    assert data == [[0, 1], [2, 3], [4]]


def test_map_chembl_ids_to_uniprot(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify rate limiter and mapping integration."""
    calls: list[int] = []

    class DummyLimiter:
        def acquire(self) -> None:
            calls.append(1)

    def fake_get_limiter(
        name: str, rps: float, burst: int | None = None
    ) -> DummyLimiter:
        return DummyLimiter()

    def fake_map(cid: str, cfg: UniprotMappingCfg) -> str:
        return f"UNI_{cid}"

    monkeypatch.setattr(mbl, "get_limiter", fake_get_limiter)
    monkeypatch.setattr(mbl, "map_chembl_to_uniprot", fake_map)

    cfg = UniprotMappingCfg()
    ids = ["CHEMBL1", "CHEMBL2", "CHEMBL3"]
    result = mbl.map_chembl_ids_to_uniprot(
        ids, cfg, batch_size=2, rps=10.0, max_workers=2
    )

    assert result == {cid: f"UNI_{cid}" for cid in ids}
    assert len(calls) == len(ids)


def test_map_chembl_ids_to_uniprot_logs_failures(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    buffer = io.StringIO()
    configure_logger(LoggerConfig(level="WARNING", stream=buffer))

    class DummyLimiter:
        def acquire(self) -> None:  # pragma: no cover - behaviour is trivial
            pass

    def fake_get_limiter(
        name: str, rps: float, burst: int | None = None
    ) -> DummyLimiter:
        return DummyLimiter()

    def failing_map(cid: str, cfg: UniprotMappingCfg) -> str:
        raise RuntimeError("boom")

    monkeypatch.setattr(mbl, "get_limiter", fake_get_limiter)
    monkeypatch.setattr(mbl, "map_chembl_to_uniprot", failing_map)

    cfg = UniprotMappingCfg()
    try:
        result = mbl.map_chembl_ids_to_uniprot(
            ["CHEMBL1"], cfg, batch_size=1, rps=1.0, max_workers=1
        )
    finally:
        configure_logger(LoggerConfig(stream=sys.stdout))

    assert result == {"CHEMBL1": None}
    records = [json.loads(line) for line in buffer.getvalue().splitlines() if line]
    failed = next(record for record in records if record["event"] == "batch_map_failed")
    assert failed["chembl_id"] == "CHEMBL1"
    assert failed["error"] == "boom"

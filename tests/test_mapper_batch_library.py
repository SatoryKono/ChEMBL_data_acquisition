from __future__ import annotations

import pytest

import library.clients.mapper_batch as mbl
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

from __future__ import annotations

from collections.abc import Generator

import pandas as pd
import pytest

from library.config import ApiCfg, TestitemBatchRetryCfg
from library.pipelines.testitem import cli


class _StreamingChemblLib:
    def __init__(self) -> None:
        self.calls: list[tuple[tuple[str, ...], dict[str, object]]] = []

    def get_testitem(self, batch: list[str], **kwargs: object) -> pd.DataFrame:
        self.calls.append((tuple(batch), kwargs))
        return pd.DataFrame(
            {"molecule_chembl_id": pd.Series(batch, dtype="string")}
        )


class _SentinelClient:
    def __repr__(self) -> str:  # pragma: no cover - debug helper
        return "<SentinelClient>"


def _large_identifier_stream(count: int) -> Generator[str, None, None]:
    for index in range(count):
        yield f"CHEMBL{index}"


@pytest.mark.slow
@pytest.mark.pipeline_scenario("degradation")
def test_fetch_testitems__peak_memory_under_control(monkeypatch: pytest.MonkeyPatch) -> None:
    psutil = pytest.importorskip("psutil")

    fake_library = _StreamingChemblLib()
    monkeypatch.setattr(cli, "_load_chembl_library", lambda: fake_library)

    identifier_count = 100_000
    ids_iter = _large_identifier_stream(identifier_count)

    status, chunk_iter, requested_unique = cli.fetch_testitems(
        ids_iter,
        api_cfg=ApiCfg(),
        batch_size=5_000,
        timeout=1.0,
        client=_SentinelClient(),
        fields=None,
        page_limit=1_000,
        retry_cfg=TestitemBatchRetryCfg(enable=False),
    )

    assert status == 0
    assert chunk_iter is not None

    process = psutil.Process()
    baseline_rss = process.memory_info().rss
    peak_rss = baseline_rss
    rows_seen = 0

    for chunk in chunk_iter:
        rows_seen += len(chunk)
        peak_rss = max(peak_rss, process.memory_info().rss)

    assert rows_seen == identifier_count
    assert len(fake_library.calls) == identifier_count // 5_000
    assert len(requested_unique) == identifier_count

    peak_increase_mb = (peak_rss - baseline_rss) / (1024 * 1024)
    assert peak_increase_mb < 128, f"Peak RSS grew by {peak_increase_mb:.2f} MiB"

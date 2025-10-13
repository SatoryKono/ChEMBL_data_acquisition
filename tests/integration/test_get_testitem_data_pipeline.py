"""Integration tests covering logging and metadata for the test item pipeline."""

from __future__ import annotations

import gc
import types
from collections.abc import Iterable

import pandas as pd
import pytest

from library.config import ApiCfg, TestitemBatchRetryCfg
from library.pipelines.testitem import cli


class _RecordingLogger:
    """Capture structured log events emitted by the pipeline module."""

    def __init__(self) -> None:
        self.events: list[tuple[str, str, dict[str, object]]] = []

    def _record(
        self,
        level: str,
        event: str,
        args: tuple[object, ...],
        kwargs: dict[str, object],
    ) -> None:
        payload = dict(kwargs)
        if args:
            payload["args"] = args
        self.events.append((level, event, payload))

    def debug(self, event: str, *args: object, **kwargs: object) -> None:
        self._record("debug", event, args, dict(kwargs))

    def info(self, event: str, *args: object, **kwargs: object) -> None:
        self._record("info", event, args, dict(kwargs))

    def warning(self, event: str, *args: object, **kwargs: object) -> None:
        self._record("warning", event, args, dict(kwargs))

    def error(self, event: str, *args: object, **kwargs: object) -> None:
        self._record("error", event, args, dict(kwargs))


class _FakeChemblLib:
    """Stand-in ChemBL client returning a predefined frame."""

    def __init__(self, frame: pd.DataFrame) -> None:
        self._frame = frame

    def get_testitem(
        self,
        _: Iterable[str],
        *,
        cfg: ApiCfg,
        client: object,
        chunk_size: int,
        timeout: float,
        fields: Iterable[str] | None,
        page_limit: int,
    ) -> pd.DataFrame:  # pragma: no cover - delegated to pandas copy
        return self._frame.copy()


class _SentinelClient:
    """Minimal ChemBL client placeholder used in tests."""

    def __repr__(self) -> str:  # pragma: no cover - debug helper
        return "<SentinelClient>"


def _run_fetch(
    monkeypatch: pytest.MonkeyPatch,
    logger: _RecordingLogger,
    frame: pd.DataFrame,
    identifiers: list[str],
) -> list[pd.DataFrame]:
    """Execute ``fetch_testitems`` with ``frame`` as the API response."""

    monkeypatch.setattr(cli, "logger", logger)
    monkeypatch.setattr(cli, "_load_chembl_library", lambda: _FakeChemblLib(frame))

    status, chunk_iter, _ = cli.fetch_testitems(
        identifiers,
        api_cfg=ApiCfg(),
        batch_size=max(len(identifiers), 1),
        timeout=1.0,
        client=_SentinelClient(),
        sample_ids=tuple(identifiers[:1]),
        fields=(),
        page_limit=50,
        retry_cfg=TestitemBatchRetryCfg(enable=False),
    )

    assert status == 0
    assert chunk_iter is not None
    return list(chunk_iter)


def _make_duplicate_frame(identifiers: Iterable[str]) -> pd.DataFrame:
    """Return a frame containing duplicate rows for each identifier."""

    duplicated = [identifier for identifier in identifiers for _ in range(2)]
    return pd.DataFrame({"molecule_chembl_id": duplicated}, dtype="string")


def _find_duplicate_event(logger: _RecordingLogger) -> dict[str, object]:
    """Return the payload of the duplicate warning emitted during fetch."""

    for level, event, payload in logger.events:
        if level == "warning" and event == "chembl_duplicate_identifiers":
            return payload
    pytest.fail("chembl_duplicate_identifiers event not emitted")


@pytest.mark.integration
@pytest.mark.pipeline_scenario("logging")
def test_fetch_testitems_memory__streaming_generator(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Large streaming inputs do not cause significant RSS growth."""

    try:
        import psutil  # type: ignore
    except ImportError:  # pragma: no cover - optional dependency fallback
        import resource

        class _Process:
            def memory_info(self) -> types.SimpleNamespace:
                usage = resource.getrusage(resource.RUSAGE_SELF)
                rss_bytes = int(usage.ru_maxrss) * 1024
                return types.SimpleNamespace(rss=rss_bytes)

        psutil = types.SimpleNamespace(Process=lambda: _Process())

    gc.collect()
    process = psutil.Process()
    baseline_rss = process.memory_info().rss

    frame = pd.DataFrame({"molecule_chembl_id": ["CHEMBL0"]}, dtype="string")
    monkeypatch.setattr(cli, "_load_chembl_library", lambda: _FakeChemblLib(frame))

    def _identifier_stream() -> Iterable[str]:
        for index in range(100_000):
            yield f"CHEMBL{index:06d}"

    sample_ids = tuple(f"CHEMBL{index:06d}" for index in range(cli._FETCH_ERROR_SAMPLE_SIZE))

    status, chunk_iter, requested = cli.fetch_testitems(
        _identifier_stream(),
        api_cfg=ApiCfg(),
        batch_size=4096,
        timeout=1.0,
        client=_SentinelClient(),
        sample_ids=sample_ids,
        fields=(),
        page_limit=50,
        retry_cfg=TestitemBatchRetryCfg(enable=False),
    )

    assert status == 0
    assert chunk_iter is not None

    peak_rss = baseline_rss
    for chunk in chunk_iter:
        assert isinstance(chunk, pd.DataFrame)
        peak_rss = max(peak_rss, process.memory_info().rss)

    final_rss = process.memory_info().rss
    assert len(requested) == 100_000

    rss_increase = max(peak_rss, final_rss) - baseline_rss
    assert rss_increase < 128 * 1024 * 1024


@pytest.mark.integration
@pytest.mark.pipeline_scenario("logging")
def test_fetch_testitems_logging__duplicate_sample_not_truncated(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Duplicate warnings expose the complete set when below the sample size."""

    logger = _RecordingLogger()
    sample_identifiers = ["CHEMBL0001", "CHEMBL0002", "CHEMBL0003"]
    frame = _make_duplicate_frame(sample_identifiers)

    chunks = _run_fetch(monkeypatch, logger, frame, sample_identifiers)
    assert len(chunks) == 1

    event = _find_duplicate_event(logger)
    assert event["duplicate_count"] == len(sample_identifiers)
    assert event["duplicate_ids"] == sample_identifiers
    assert event["duplicates_truncated"] is False


@pytest.mark.integration
@pytest.mark.pipeline_scenario("logging")
def test_fetch_testitems_logging__duplicate_sample_truncated(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Duplicate warnings include truncated samples when exceeding the limit."""

    logger = _RecordingLogger()
    sample_size = cli._DUPLICATE_IDENTIFIER_LOG_SAMPLE_SIZE
    identifiers = [
        f"CHEMBL{index:04d}"
        for index in range(sample_size + 3)
    ]
    frame = _make_duplicate_frame(identifiers)

    chunks = _run_fetch(monkeypatch, logger, frame, identifiers)
    assert len(chunks) == 1

    event = _find_duplicate_event(logger)
    assert event["duplicate_count"] == len(identifiers)
    assert len(event["duplicate_ids"]) == sample_size
    assert event["duplicate_ids"] == identifiers[:sample_size]
    assert event["duplicates_truncated"] is True

from __future__ import annotations

from pathlib import Path
from typing import Any
import hashlib

import pandas as pd
import pytest
import yaml

from library import pubmed_library as pl
from library.clients import pubmed as pubmed_client
from library.config import RetryCfg


def test_pubmed_library_main_smoke(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Run the CLI with stubbed network calls."""
    input_csv = Path("tests/data/pmids.csv")
    output_csv = tmp_path / "out.csv"
    verbose_output_csv = tmp_path / "out_verbose.csv"

    seen_retry_cfg: list[RetryCfg | None] = []

    def fake_fetch_pubmed_batch(
        session,
        pmids,
        delay,
        cfg=None,
        *,
        retry_cfg: RetryCfg | None = None,
        client=None,
    ):
        seen_retry_cfg.append(retry_cfg)
        return [
            {
                "PubMed.PMID": pid if pid != "2" else "",
                "PubMed.DOI": "10.1000/example",
                "PubMed.ArticleTitle": "Example",
            }
            for pid in pmids
        ]

    def fake_fetch_semantic_scholar_batch(session, pmids, delay, cfg=None):
        return [
            {
                "scholar.PMID": pid,
                "scholar.DOI": "10.1000/example",
                "scholar.Error": "",
            }
            for pid in pmids
        ]

    def fake_fetch_openalex(session, pmid, cfg=None, limiter=None):
        return {"OpenAlex.Error": ""}

    def fake_fetch_crossref(session, doi, cfg=None, limiter=None):
        return {"crossref.Error": ""}

    class DummyLimiter:
        def acquire(self) -> None:
            return None

    monkeypatch.setattr(pl, "fetch_pubmed_batch", fake_fetch_pubmed_batch)
    monkeypatch.setattr(
        pl, "fetch_semantic_scholar_batch", fake_fetch_semantic_scholar_batch
    )
    monkeypatch.setattr(pl, "fetch_openalex", fake_fetch_openalex)
    monkeypatch.setattr(pl, "fetch_crossref", fake_fetch_crossref)
    monkeypatch.setattr(pl, "get_limiter", lambda *args, **kwargs: DummyLimiter())

    levels: list[str] = []

    def fake_print_results(records, *, level: str = "DEBUG") -> None:
        levels.extend([level] * len(records))

    monkeypatch.setattr(pl, "print_results", fake_print_results)

    exit_code = pl.main(
        [
            "--input-csv",
            str(input_csv),
            "--output",
            str(output_csv),
            "--log-level",
            "INFO",
        ]
    )
    assert exit_code == 1
    assert levels
    assert all(level == "DEBUG" for level in levels)
    assert output_csv.exists()
    df = pd.read_csv(output_csv)
    assert not df.empty
    failure_path = output_csv.with_name(f"{output_csv.stem}_failure_cases.csv")
    assert failure_path.exists()
    failure_df = pd.read_csv(failure_path)
    assert not failure_df.empty

    quality_base = output_csv.with_suffix("")
    quality_report = quality_base.with_name(
        f"{quality_base.name}_quality_report_table.csv"
    )
    correlation_report = quality_base.with_name(
        f"{quality_base.name}_data_correlation_report_table.csv"
    )
    assert quality_report.exists()
    assert correlation_report.exists()

    meta_path = output_csv.with_name(output_csv.name + ".meta.yaml")
    assert meta_path.exists()
    digest = hashlib.sha256(output_csv.read_bytes()).hexdigest()
    meta = yaml.safe_load(meta_path.read_text())
    assert meta["stats"]["output_sha256"] == digest
    assert meta["schema"] == "PubMedDocumentsSchema"

    start_idx = len(levels)
    second_output = tmp_path / "out_second.csv"
    exit_code_verbose = pl.main(
        [
            "--input-csv",
            str(input_csv),
            "--output",
            str(verbose_output_csv),
            "--log-level",
            "INFO",
            "--keep-verbose-dumps",
        ]
    )
    assert exit_code_verbose == 1
    verbose_levels = levels[start_idx:]
    assert verbose_levels
    assert all(level == "INFO" for level in verbose_levels)

    second_exit_code = pl.main(
        [
            "--input-csv",
            str(input_csv),
            "--output",
            str(second_output),
            "--log-level",
            "INFO",
        ]
    )
    assert second_exit_code == 1
    assert second_output.exists()
    assert output_csv.read_bytes() == second_output.read_bytes()
    second_digest = hashlib.sha256(second_output.read_bytes()).hexdigest()
    assert second_digest == digest

    second_meta_path = second_output.with_name(second_output.name + ".meta.yaml")
    assert second_meta_path.exists()
    second_meta = yaml.safe_load(second_meta_path.read_text())
    assert second_meta["stats"]["output_sha256"] == second_digest
    assert second_meta["schema"] == "PubMedDocumentsSchema"

    assert seen_retry_cfg
    assert all(cfg is not None for cfg in seen_retry_cfg)
    expected_backoff = pl.Config().retry.backoff_factor
    assert all(
        cfg.backoff_factor == pytest.approx(expected_backoff)
        for cfg in seen_retry_cfg
        if cfg is not None
    )


def test_pubmed_client_retry_logging(monkeypatch: pytest.MonkeyPatch) -> None:
    """PubMed client should log retry delay with jitter for 429 responses."""

    responses = iter(
        [
            (429, "Too Many Requests", None, "", {"Retry-After": "1.7"}),
            (200, "<xml/>", "<xml/>", "", {}),
        ]
    )

    def fake_make_request(
        *_args: Any, **_kwargs: Any
    ) -> tuple[int, str, Any, str, dict[str, str]]:
        try:
            return next(responses)
        except StopIteration:  # pragma: no cover - defensive
            raise AssertionError("Unexpected request count")

    sleep_calls: list[float] = []

    class RecordingLogger:
        def __init__(self) -> None:
            self.records: list[tuple[str, str, dict[str, Any]]] = []

        def info(
            self, event: str, *args: Any, extra: dict[str, Any] | None = None, **kwargs: Any
        ) -> None:
            self.records.append(("info", event, extra.copy() if extra else {}))

        def debug(
            self, event: str, *args: Any, extra: dict[str, Any] | None = None, **kwargs: Any
        ) -> None:
            self.records.append(("debug", event, extra.copy() if extra else {}))

    monkeypatch.setattr(pubmed_client, "_make_request", fake_make_request)
    monkeypatch.setattr(pubmed_client, "sleep", lambda value: sleep_calls.append(value))
    monkeypatch.setattr(pubmed_client.random, "uniform", lambda _a, _b: 0.2)

    fake_logger = RecordingLogger()
    monkeypatch.setattr(pubmed_client, "logger", fake_logger)

    retry_cfg = RetryCfg(max_attempts=2, backoff_factor=0.5)
    data, error = pubmed_client._do_request(
        session=object(),
        url="http://example",
        delay=0.3,
        expect_json=False,
        retries=1,
        timeout=(1, 5),
        retry_cfg=retry_cfg,
    )

    assert data == "<xml/>"
    assert error == ""

    base = retry_cfg.backoff_factor * (2 ** (1 - 1))
    expected_delay = pytest.approx(1.7)
    assert sleep_calls == [expected_delay]

    retry_records = [r for r in fake_logger.records if r[1] == "request_retry"]
    assert retry_records and "delay" in retry_records[0][2]
    assert retry_records[0][2]["delay"] == expected_delay

    sleep_logs = [r for r in fake_logger.records if r[1] == "retry_sleep"]
    assert sleep_logs and sleep_logs[0][2]["delay"] == expected_delay

    responses_no_header = iter(
        [
            (429, "Too Many Requests", None, "", {}),
            (200, "<xml/>", "<xml/>", "", {}),
        ]
    )

    def fake_make_request_no_retry_after(
        *_args: Any, **_kwargs: Any
    ) -> tuple[int, str, Any, str, dict[str, str]]:
        try:
            return next(responses_no_header)
        except StopIteration:  # pragma: no cover - defensive
            raise AssertionError("Unexpected request count for fallback scenario")

    monkeypatch.setattr(pubmed_client, "_make_request", fake_make_request_no_retry_after)
    jitter_value = 0.3
    monkeypatch.setattr(pubmed_client.random, "uniform", lambda _a, _b: jitter_value)

    sleep_calls.clear()
    fallback_logger = RecordingLogger()
    monkeypatch.setattr(pubmed_client, "logger", fallback_logger)

    base_delay = 0.0
    data_retry, error_retry = pubmed_client._do_request(
        session=object(),
        url="http://example",
        delay=base_delay,
        expect_json=False,
        retries=1,
        timeout=(1, 5),
        retry_cfg=None,
    )

    assert data_retry == "<xml/>"
    assert error_retry == ""

    fallback_retry_records = [
        r for r in fallback_logger.records if r[1] == "request_retry"
    ]
    assert fallback_retry_records and "delay" in fallback_retry_records[0][2]
    fallback_delay = fallback_retry_records[0][2]["delay"]
    default_backoff = RetryCfg().backoff_factor * (2 ** (1 - 1))
    timeout_cap = pubmed_client._max_timeout((1, 5))
    assert timeout_cap is not None
    expected_fallback_delay = min(
        max(base_delay, default_backoff + jitter_value),
        timeout_cap,
    )
    assert fallback_delay == pytest.approx(expected_fallback_delay)
    assert fallback_delay > 0

    fallback_sleep_logs = [
        r for r in fallback_logger.records if r[1] == "retry_sleep"
    ]
    assert fallback_sleep_logs and fallback_sleep_logs[0][2]["delay"] == fallback_delay
    assert sleep_calls == [pytest.approx(fallback_delay)]
    assert sleep_calls[0] > 0

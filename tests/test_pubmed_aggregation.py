"""Tests for :mod:`library.pubmed.aggregation`."""

from __future__ import annotations

import logging

import pytest

from library.pubmed import aggregation as pa


def test_merge_records() -> None:
    a = {"a": 1, "b": 2}
    b = {"b": 3, "c": 4}
    merged = pa.merge_records(a, b)
    assert merged == {"a": 1, "b": 3, "c": 4}


def test_print_results_logs_at_debug_by_default(
    caplog: pytest.LogCaptureFixture, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = logging.getLogger("test.debug")
    monkeypatch.setattr(pa, "logger", fake)
    rec = {"PubMed.ArticleTitle": "Example title"}

    with caplog.at_level(logging.INFO, logger="test.debug"):
        pa.print_results([rec])
    assert all("Example title" not in record.message for record in caplog.records)

    caplog.clear()
    with caplog.at_level(logging.DEBUG, logger="test.debug"):
        pa.print_results([rec])
    assert any("Example title" in record.message for record in caplog.records)


def test_print_results_can_log_at_info(
    caplog: pytest.LogCaptureFixture, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = logging.getLogger("test.info")
    monkeypatch.setattr(pa, "logger", fake)
    rec = {"PubMed.ArticleTitle": "Example title"}

    with caplog.at_level(logging.INFO, logger="test.info"):
        pa.print_results([rec], level="INFO")

    assert any("Example title" in record.message for record in caplog.records)

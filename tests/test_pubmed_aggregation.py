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


def test_print_results(
    caplog: pytest.LogCaptureFixture, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = logging.getLogger("test")
    caplog.set_level("INFO", logger="test")
    monkeypatch.setattr(pa, "logger", fake)
    rec = {"PubMed.ArticleTitle": "Example title"}
    pa.print_results([rec])
    assert "Example title" in caplog.text

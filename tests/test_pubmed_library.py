"""Tests for :mod:`library.pubmed_library`."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import requests
import pytest

from library import pubmed_library as pl
from library.config import CrossRefCfg, OpenAlexCfg


DATA_DIR = Path(__file__).parent / "data"


def test_read_pmids() -> None:
    path = DATA_DIR / "pmids.csv"
    pmids = pl.read_pmids(path)
    assert pmids == ["1", "2"]


def test_fetch_openalex_uses_cfg(monkeypatch) -> None:
    cfg = OpenAlexCfg(
        base="https://example.org",
        timeout_connect=1,
        timeout_read=2,
        retries=4,
        rps=10,
        mailto="info@example.org",
    )

    captured: Dict[str, Any] = {}

    def fake_do_request(
        session,
        url,
        sleep,
        expect_json=True,
        retries=2,
        method="GET",
        timeout=10,
        **kwargs,
    ):
        captured.update(
            {"url": url, "sleep": sleep, "retries": retries, "timeout": timeout}
        )
        return {
            "type": "article",
            "type_crossref": "journal-article",
            "genre": "journal",
            "id": "123",
            "host_venue": {"display_name": "Venue"},
            "mesh": [],
        }, ""

    sleeps: list[float] = []
    monkeypatch.setattr(pl, "_do_request", fake_do_request)
    monkeypatch.setattr(pl.time, "sleep", lambda s: sleeps.append(s))

    session = requests.Session()
    res = pl.fetch_openalex(session, "1", cfg=cfg)
    assert res["OpenAlex.Error"] == ""
    assert (
        captured["url"] == "https://example.org/works/pmid:1?mailto=info%40example.org"
    )
    assert captured["timeout"] == (1, 2)
    assert captured["retries"] == 4
    assert captured["sleep"] == pytest.approx(0.1)
    assert sleeps == [pytest.approx(0.1)]


def test_fetch_crossref_uses_cfg(monkeypatch) -> None:
    cfg = CrossRefCfg(
        base="https://api.example.org",
        timeout_connect=2,
        timeout_read=3,
        retries=5,
        rps=8,
        mailto="info@example.org",
    )

    captured: Dict[str, Any] = {}

    def fake_do_request(
        session,
        url,
        sleep,
        expect_json=True,
        retries=2,
        method="GET",
        timeout=10,
        **kwargs,
    ):
        captured.update(
            {"url": url, "sleep": sleep, "retries": retries, "timeout": timeout}
        )
        return {"message": {}}, ""

    sleeps: list[float] = []
    monkeypatch.setattr(pl, "_do_request", fake_do_request)
    monkeypatch.setattr(pl.time, "sleep", lambda s: sleeps.append(s))

    session = requests.Session()
    res = pl.fetch_crossref(session, "10.1000/xyz", cfg=cfg)
    assert res["crossref.Error"] == ""
    assert (
        captured["url"]
        == "https://api.example.org/works/10.1000%2Fxyz?mailto=info%40example.org"
    )
    assert captured["timeout"] == (2, 3)
    assert captured["retries"] == 5
    assert captured["sleep"] == pytest.approx(1 / 8)
    assert sleeps == [pytest.approx(1 / 8)]

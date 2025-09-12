"""Tests for :mod:`library.uniprot_library`."""

from __future__ import annotations

import time

import pytest

requests = pytest.importorskip("requests")
responses = pytest.importorskip("responses")

from library import uniprot_library as ul  # noqa: E402
from library.config import UniprotCfg  # noqa: E402


def test_extract_names() -> None:
    sample = {
        "proteinDescription": {
            "recommendedName": {"fullName": {"value": "Protein X"}},
            "alternativeNames": [{"fullName": {"value": "Alt Name"}}],
        },
        "genes": [{"geneName": {"value": "GENE1"}, "synonyms": [{"value": "G1"}]}],
    }
    names = ul.extract_names(sample)
    assert names == {"Protein X", "Alt Name", "GENE1", "G1"}


@responses.activate
def test_fetch_uniprot_network_error() -> None:
    cfg = UniprotCfg(base="https://example.org", delay=0)
    responses.add(
        responses.GET,
        "https://example.org/uniprotkb/P12345.json",
        body=requests.RequestException("boom"),
    )
    with pytest.raises(ul.UniProtFetchError):
        ul.fetch_uniprot("P12345", cfg=cfg)


@responses.activate
def test_fetch_uniprot_bad_json() -> None:
    cfg = UniprotCfg(base="https://example.org", delay=0)
    responses.add(
        responses.GET,
        "https://example.org/uniprotkb/P12345.json",
        body="{",
        status=200,
    )
    with pytest.raises(ul.UniProtFetchError):
        ul.fetch_uniprot("P12345", cfg=cfg)


@responses.activate
def test_fetch_uniprot_uses_cfg(monkeypatch) -> None:
    called: dict[str, object] = {}
    orig_get = ul._session.get

    def capture(url: str, timeout: tuple[int, int]):
        called["url"] = url
        called["timeout"] = timeout
        return orig_get(url, timeout=timeout)

    sleeps: list[float] = []

    monkeypatch.setattr(ul._session, "get", capture)
    monkeypatch.setattr(time, "sleep", lambda s: sleeps.append(s))
    cfg = UniprotCfg(
        base="https://example.org/api",
        timeout_connect=1,
        timeout_read=2,
        delay=0.5,
    )
    responses.add(
        responses.GET, "https://example.org/api/uniprotkb/P12345.json", json={}
    )
    ul.fetch_uniprot("P12345", cfg=cfg)
    assert called["url"] == "https://example.org/api/uniprotkb/P12345.json"
    assert called["timeout"] == (1, 2)
    assert sleeps and sleeps[0] == pytest.approx(0.5)

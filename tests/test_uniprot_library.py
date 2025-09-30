"""Tests for :mod:`library.uniprot_library`."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

requests = pytest.importorskip("requests")
responses = pytest.importorskip("responses")

from library import uniprot_library as ul  # noqa: E402
from library.clients import uniprot as uniprot_client  # noqa: E402
from library.config import IupharCfg, RetryCfg, UniprotCfg  # noqa: E402


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
def test_fetch_uniprot_network_error(monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = UniprotCfg(base="https://example.org", delay=0)
    retry_cfg = RetryCfg(max_attempts=2, backoff_factor=0.1)
    monkeypatch.setattr(uniprot_client, "_retry_cfg", retry_cfg)
    sleeps: list[float] = []
    monkeypatch.setattr(uniprot_client, "sleep", lambda seconds: sleeps.append(seconds))
    monkeypatch.setattr(uniprot_client.random, "uniform", lambda _a, _b: 0.0)
    responses.add(
        responses.GET,
        "https://example.org/uniprotkb/P12345.json",
        body=requests.RequestException("boom"),
        repeat=True,
    )
    with pytest.raises(ul.UniProtFetchError):
        ul.fetch_uniprot("P12345", cfg=cfg)
    assert sleeps == [pytest.approx(0.1)]


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
def test_fetch_uniprot_uses_cfg(monkeypatch: pytest.MonkeyPatch) -> None:
    called: dict[str, object] = {}
    session = uniprot_client.get_session()
    orig_get = session.get

    def capture(url: str, timeout: tuple[int, int]):
        called["url"] = url
        called["timeout"] = timeout
        return orig_get(url, timeout=timeout)

    monkeypatch.setattr(session, "get", capture)
    sleeps: list[float] = []
    monkeypatch.setattr(uniprot_client, "sleep", lambda seconds: sleeps.append(seconds))
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
    assert sleeps == []


@responses.activate
def test_fetch_uniprot_retries_with_backoff(monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = UniprotCfg(base="https://example.org", delay=0.5)
    retry_cfg = RetryCfg(max_attempts=3, backoff_factor=0.2)
    monkeypatch.setattr(uniprot_client, "_retry_cfg", retry_cfg)
    limiter_calls: list[tuple[int, int]] = []
    acquire_calls = 0

    class DummyLimiter:
        def acquire(self) -> None:
            nonlocal acquire_calls
            acquire_calls += 1

    def fake_get_limiter(name: str, rps: int, burst: int) -> DummyLimiter:
        limiter_calls.append((rps, burst))
        return DummyLimiter()

    monkeypatch.setattr(uniprot_client, "get_limiter", fake_get_limiter)
    sleeps: list[float] = []
    monkeypatch.setattr(uniprot_client, "sleep", lambda seconds: sleeps.append(seconds))
    monkeypatch.setattr(uniprot_client.random, "uniform", lambda _a, _b: 0.0)

    responses.add(
        responses.GET,
        "https://example.org/uniprotkb/P12345.json",
        status=500,
    )
    responses.add(
        responses.GET,
        "https://example.org/uniprotkb/P12345.json",
        json={"primaryAccession": "P12345"},
    )

    result = ul.fetch_uniprot("P12345", cfg=cfg)

    assert result == {"primaryAccession": "P12345"}
    assert limiter_calls == [(cfg.rps, cfg.burst), (cfg.rps, cfg.burst)]
    assert acquire_calls == 2
    assert sleeps == [pytest.approx(0.7)]


@responses.activate
def test_collect_info_enriches_gtop(tmp_path: Path) -> None:
    data = {
        "uniProtKBCrossReferences": [
            {"database": "GuidetoPHARMACOLOGY", "id": "1234"},
            {"database": "GuidetoPHARMACOLOGY", "id": "5678"},
        ]
    }
    data_dir = tmp_path
    (data_dir / "P12345.json").write_text(json.dumps(data))

    gtop_cfg = IupharCfg(base="https://gtop.example.org/services", rps=10, burst=10)
    cfg = UniprotCfg(base="https://example.org", delay=0)

    responses.add(
        responses.GET,
        "https://gtop.example.org/services/targets/1234/naturalLigands",
        json=[{"ligandId": 1}, {"ligandId": 2}],
    )
    responses.add(
        responses.GET,
        "https://gtop.example.org/services/targets/1234/interactions",
        json=[{"interactionId": 1}, {"interactionId": 2}, {"interactionId": 3}],
    )
    responses.add(
        responses.GET,
        "https://gtop.example.org/services/targets/1234/function",
        json=[
            {
                "description": "Physiological function",
                "property": "Regulates sample process",
            }
        ],
    )

    result = ul.collect_info("P12345", data_dir=data_dir, cfg=cfg, gtop_cfg=gtop_cfg)

    assert result["gtop_natural_ligands_n"] == "2"
    assert result["gtop_interactions_n"] == "3"
    assert (
        result["gtop_function_text_short"]
        == "Physiological function: Regulates sample process"
    )
    assert len(responses.calls) == 3

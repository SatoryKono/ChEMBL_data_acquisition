"""Tests for :mod:`library.uniprot_library`."""

from __future__ import annotations

import json
import time
from pathlib import Path

import pytest

requests = pytest.importorskip("requests")
responses = pytest.importorskip("responses")

from library import uniprot_library as ul  # noqa: E402
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
def test_fetch_uniprot_network_error(monkeypatch) -> None:
    cfg = UniprotCfg(base="https://example.org", delay=0)
    monkeypatch.setattr(ul, "_retry_cfg", RetryCfg(max_attempts=1))
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


@responses.activate
def test_fetch_uniprot_retries_transient_error(monkeypatch) -> None:
    cfg = UniprotCfg(base="https://example.org", delay=0)
    monkeypatch.setattr(
        ul,
        "_retry_cfg",
        RetryCfg(max_attempts=3, backoff_factor=0.5, status_forcelist=[500]),
    )
    url = "https://example.org/uniprotkb/P12345.json"

    responses.add(responses.GET, url, status=500)
    responses.add(responses.GET, url, json={"status": "ok"})

    sleeps: list[float] = []
    monkeypatch.setattr(ul, "sleep", lambda delay: sleeps.append(delay))
    monkeypatch.setattr(ul.random, "uniform", lambda _a, b: b)

    records: list[tuple[str, dict[str, object]]] = []

    def capture(event: str, **payload: object) -> None:
        records.append((event, payload))

    monkeypatch.setattr(ul.logger, "warning", capture)

    result = ul.fetch_uniprot("P12345", cfg=cfg)

    assert result == {"status": "ok"}
    assert len(responses.calls) == 2
    assert sleeps == [pytest.approx(1.0)]
    assert len(records) == 1
    event, payload = records[0]
    assert event == "uniprot_request_retry"
    assert payload["url"] == url
    assert payload["attempt"] == 1
    assert payload["max_attempts"] == 3
    assert payload["status"] == 500
    assert payload["rps"] == cfg.rps
    assert "500" in str(payload["error"])


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

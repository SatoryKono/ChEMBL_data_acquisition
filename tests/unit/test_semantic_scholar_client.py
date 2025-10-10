from __future__ import annotations

import pytest
import requests

from library.clients import semantic_scholar
from library.config.models import SemanticScholarCfg


@pytest.mark.unit
def test_fetch_semantic_scholar__includes_api_key_header(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def _fake_do_request(session, url, delay, *, headers, **kwargs):  # type: ignore[override]
        captured["headers"] = headers
        return (
            {
                "paperId": "S2:123",
                "externalIds": {"DOI": "10.1234/example"},
                "publicationTypes": ["Journal Article"],
                "venue": "Example Journal",
            },
            "",
        )

    monkeypatch.setattr(semantic_scholar, "_do_request", _fake_do_request)

    session = requests.Session()
    cfg = SemanticScholarCfg(api_key="  secret-token  ")

    result = semantic_scholar.fetch_semantic_scholar(session, "123456", 0.0, cfg=cfg)

    assert result["scholar.Error"] == ""
    headers = captured.get("headers")
    assert isinstance(headers, dict)
    assert headers["Accept"] == "application/json"
    assert headers["x-api-key"] == "secret-token"


@pytest.mark.unit
def test_fetch_semantic_scholar_batch__omits_api_key_header_when_missing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def _fake_do_request(session, url, delay, *, headers, **kwargs):  # type: ignore[override]
        captured["headers"] = headers
        return (
            [
                {
                    "paperId": "S2:123",
                    "externalIds": {"DOI": "10.1234/example"},
                    "publicationTypes": ["Journal Article"],
                    "venue": "Example Journal",
                    "pmid": "123456",
                }
            ],
            "",
        )

    monkeypatch.setattr(semantic_scholar, "_do_request", _fake_do_request)

    session = requests.Session()

    results = semantic_scholar.fetch_semantic_scholar_batch(session, ["123456"], 0.0)

    assert results
    assert results[0]["scholar.Error"] == ""
    headers = captured.get("headers")
    assert isinstance(headers, dict)
    assert headers["Accept"] == "application/json"
    assert "x-api-key" not in headers


@pytest.mark.unit
def test_semantic_scholar_cfg__rejects_blank_api_key() -> None:
    with pytest.raises(ValueError):
        SemanticScholarCfg(api_key="   ")

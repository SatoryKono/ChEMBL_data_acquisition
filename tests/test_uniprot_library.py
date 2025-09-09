"""Tests for :mod:`library.uniprot_library`."""

from __future__ import annotations

import json
import requests
import pytest

from library import uniprot_library as ul


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


class _DummyResponse:
    def __init__(self, *, json_exc: Exception | None = None) -> None:
        self.status_code = 200
        self._json_exc = json_exc
        self.text = "{}"

    def raise_for_status(self) -> None:  # pragma: no cover - always ok
        return None

    def json(self) -> dict[str, str]:
        if self._json_exc:
            raise self._json_exc
        return {}

    def close(self) -> None:
        return None

    def __enter__(self) -> "_DummyResponse":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # pragma: no cover - no cleanup
        return None


def test_fetch_uniprot_network_error(monkeypatch) -> None:
    def fake_get(url: str, timeout: float) -> None:
        raise requests.RequestException("boom")

    monkeypatch.setattr(ul._session, "get", fake_get)
    with pytest.raises(ul.UniProtFetchError):
        ul.fetch_uniprot("P12345")


def test_fetch_uniprot_bad_json(monkeypatch) -> None:
    resp = _DummyResponse(json_exc=json.JSONDecodeError("msg", "doc", 0))

    def fake_get(url: str, timeout: float) -> _DummyResponse:
        return resp

    monkeypatch.setattr(ul._session, "get", fake_get)
    with pytest.raises(ul.UniProtFetchError):
        ul.fetch_uniprot("P12345")

"""Unit tests for the ChEMBL document fetch helpers."""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd
import pytest
import requests

from library.config import ApiCfg
from library.pipelines.document.chembl_document import get_documents


class _StubChemblClient:
    """Deterministic stub that mimics :class:`ChemblClient`."""

    def __init__(self, responses: dict[str, object]) -> None:
        self._responses = responses
        self.calls: list[tuple[str, float | None]] = []

    def request_json(
        self, url: str, *, cfg: ApiCfg, timeout: float | None = None
    ) -> dict:
        del cfg  # Unused in the stub.
        self.calls.append((url, timeout))
        response = self._responses[url]
        if isinstance(response, Exception):
            raise response
        assert isinstance(response, dict)
        return response


def _chunk_url(cfg: ApiCfg, ids: Sequence[str]) -> str:
    base = cfg.chembl_base.rstrip("/")
    joined_ids = ",".join(ids)
    return f"{base}/document.json?format=json&document_chembl_id__in={joined_ids}"


def _build_response(identifier: str, title: str) -> dict[str, object]:
    return {
        "documents": [
            {
                "document_chembl_id": identifier,
                "title": title,
                "abstract": "",
                "doi": "",
                "year": 2024,
                "journal_full_title": "Journal",
                "journal": "J",
                "volume": "1",
                "issue": "2",
                "first_page": "10",
                "last_page": "20",
                "pubmed_id": "123",
                "authors": ["Author"],
            }
        ]
    }


def _build_chunk_response(identifiers: Sequence[str]) -> dict[str, object]:
    return {
        "documents": [
            {
                "document_chembl_id": identifier,
                "title": f"Title {identifier}",
                "abstract": "",
                "doi": "",
                "year": 2024,
                "journal_full_title": "Journal",
                "journal": "J",
                "volume": "1",
                "issue": "1",
                "first_page": "1",
                "last_page": "2",
                "pubmed_id": "123",
                "authors": ["Author"],
            }
            for identifier in identifiers
        ]
    }


@pytest.mark.unit
def test_get_documents__splits_chunk_on_timeout(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """The helper should split large chunks when a read timeout occurs."""

    cfg = ApiCfg(chembl_base="https://example.test/api", timeout_read=12.0)
    timeout = 9.5

    combined_url = _chunk_url(cfg, ["CHEMBL1", "CHEMBL2"])
    responses = {
        combined_url: requests.ReadTimeout("simulated timeout"),
        _chunk_url(cfg, ["CHEMBL1"]): _build_response("CHEMBL1", "Alpha"),
        _chunk_url(cfg, ["CHEMBL2"]): _build_response("CHEMBL2", "Beta"),
    }
    client = _StubChemblClient(responses)

    with caplog.at_level("WARNING"):
        frame = get_documents(
            ["CHEMBL1", "CHEMBL2"],
            cfg=cfg,
            client=client,
            chunk_size=2,
            timeout=timeout,
        )

    assert [call[0] for call in client.calls] == [
        combined_url,
        _chunk_url(cfg, ["CHEMBL1"]),
        _chunk_url(cfg, ["CHEMBL2"]),
    ]
    assert all(call[1] == timeout for call in client.calls)

    assert isinstance(frame, pd.DataFrame)
    assert list(frame["document_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]
    assert list(frame["title"]) == ["Alpha", "Beta"]

    assert any(
        record.getMessage().startswith("chembl_document_timeout_split")
        for record in caplog.records
    )


@pytest.mark.unit
def test_get_documents__propagates_timeout_for_single_identifier() -> None:
    """When the chunk is a single identifier the timeout error should bubble up."""

    cfg = ApiCfg(chembl_base="https://example.test/api", timeout_read=6.0)
    timeout = 5.0

    single_url = _chunk_url(cfg, ["CHEMBL1"])
    responses = {single_url: requests.ReadTimeout("simulated timeout")}
    client = _StubChemblClient(responses)

    with pytest.raises(requests.ReadTimeout):
        get_documents(
            ["CHEMBL1"],
            cfg=cfg,
            client=client,
            chunk_size=1,
            timeout=timeout,
        )

    assert client.calls == [(single_url, timeout)]


@pytest.mark.unit
def test_get_documents__caps_large_chunk_size(caplog: pytest.LogCaptureFixture) -> None:
    cfg = ApiCfg(chembl_base="https://example.test/api", timeout_read=60.0)
    identifiers = [f"CHEMBL{i}" for i in range(1, 27)]

    first_chunk = identifiers[:20]
    second_chunk = identifiers[20:]
    responses = {
        _chunk_url(cfg, first_chunk): _build_chunk_response(first_chunk),
        _chunk_url(cfg, second_chunk): _build_chunk_response(second_chunk),
    }
    client = _StubChemblClient(responses)

    with caplog.at_level("INFO"):
        frame = get_documents(
            identifiers,
            cfg=cfg,
            client=client,
            chunk_size=50,
            timeout=None,
        )

    assert [call[0] for call in client.calls] == [
        _chunk_url(cfg, first_chunk),
        _chunk_url(cfg, second_chunk),
    ]
    assert all(call[1] == cfg.timeout_read for call in client.calls)
    assert sorted(frame["document_chembl_id"]) == sorted(identifiers)
    assert any(
        record.getMessage().startswith("chembl_document_chunk_size_capped")
        for record in caplog.records
    )

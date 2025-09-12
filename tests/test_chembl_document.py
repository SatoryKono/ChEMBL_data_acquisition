"""Tests for :mod:`library.chembl_document`."""

from __future__ import annotations

from library.chembl_document import get_documents
from library.config import ApiCfg


class FakeClient:
    """Minimal stand-in for :class:`ChemblClient` used in tests."""

    def __init__(self) -> None:
        self.called: list[tuple[str, float | None]] = []

    def request_json(
        self, url: str, *, cfg: ApiCfg, timeout: float | None = None
    ) -> dict[str, object]:
        self.called.append((url, timeout))
        ids = url.split("document_chembl_id__in=")[-1].split("&")[0].split(",")
        docs = [
            {
                "document_chembl_id": i,
                "title": f"title-{i}",
                "abstract": "",
                "doi": "",
                "year": 2020,
                "journal": "J",
                "journal_full_title": "Journal",
                "volume": "1",
                "issue": "2",
                "first_page": "3",
                "last_page": "4",
                "pubmed_id": 1,
                "authors": "A",
            }
            for i in ids
        ]
        return {"documents": docs}


def test_get_documents_returns_frame() -> None:
    cfg = ApiCfg(
        chembl_base="https://example.com",
        user_agent="test/0.1 (mailto:test@example.com)",
    )
    client = FakeClient()
    df = get_documents(
        ["CHEMBL1", "CHEMBL2"], cfg=cfg, client=client, chunk_size=1, timeout=5.0
    )
    assert df["document_chembl_id"].tolist() == ["CHEMBL1", "CHEMBL2"]
    assert len(client.called) == 2
    assert all(t == 5.0 for _, t in client.called)


def test_get_documents_handles_invalid_ids() -> None:
    cfg = ApiCfg(
        chembl_base="https://example.com",
        user_agent="test/0.1 (mailto:test@example.com)",
    )
    client = FakeClient()
    df = get_documents(["", "#N/A"], cfg=cfg, client=client)
    assert df.empty
    assert not client.called

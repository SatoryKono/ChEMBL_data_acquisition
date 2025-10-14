"""Regression tests ensuring TLS-aware sessions are used for document lookups."""

from __future__ import annotations

import csv
from collections.abc import Callable
from pathlib import Path

import pandas as pd
import pytest
import requests

from library.config import Config
from library.pipelines.document.service import DocumentPipeline


class _StubSession(requests.Session):
    """Session subclass tracking close calls for test assertions."""

    def __init__(self, service: str) -> None:
        super().__init__()
        self.service = service
        self.closed = False
        self.headers["X-Test-Service"] = service

    def close(self) -> None:  # pragma: no cover - exercised indirectly
        self.closed = True
        super().close()


@pytest.fixture()
def _patch_document_dependencies(
    monkeypatch: pytest.MonkeyPatch,
) -> dict[str, list[_StubSession]]:
    """Stub external dependencies to isolate session orchestration."""

    calls: dict[str, list[_StubSession]] = {"openalex": [], "crossref": []}

    def _session_factory(name: str) -> Callable[..., requests.Session]:
        def _factory(*_: object, **__: object) -> requests.Session:
            session = _StubSession(name)
            calls[name].append(session)
            return session

        return _factory

    monkeypatch.setattr(
        "library.pipelines.document.service.openalex_session",
        _session_factory("openalex"),
    )
    monkeypatch.setattr(
        "library.pipelines.document.service.crossref_session",
        _session_factory("crossref"),
    )

    def _pubmed_batch(
        _session: requests.Session,
        pmids: list[str],
        _sleep: float,
        *,
        cfg: object,
        retry_cfg: object,
    ) -> list[dict[str, str]]:
        pmid = pmids[0]
        return [
            {
                "PubMed.PMID": pmid,
                "PubMed.DOI": "10.1000/example",
                "PubMed.Title": "Example title",
                "PubMed.PublicationType": "Journal Article",
                "PubMed.Error": "",
            }
        ]

    monkeypatch.setattr(
        "library.integration.pubmed_library.fetch_pubmed_batch",
        _pubmed_batch,
    )
    monkeypatch.setattr(
        "library.integration.semantic_scholar_library.fetch_semantic_scholar_batch",
        lambda *_args, **_kwargs: [],
    )

    def _semantic_single(
        _session: requests.Session,
        pmid: str,
        *,
        cfg: object,
        limiter: object | None,
        retry_cfg: object,
    ) -> dict[str, str]:
        return {
            "scholar.PMID": pmid,
            "scholar.DOI": "",
            "scholar.Error": "",
            "scholar.PublicationTypes": "",
        }

    monkeypatch.setattr(
        "library.integration.semantic_scholar_library.fetch_semantic_scholar",
        _semantic_single,
    )

    def _openalex(
        session: requests.Session,
        pmid: str,
        cfg: object,
        *,
        limiter: object | None = None,
        retry_cfg: object | None = None,
    ) -> dict[str, str]:
        assert isinstance(session, _StubSession)
        assert session.service == "openalex"
        return {
            "OpenAlex.PublicationTypes": "journal-article",
            "OpenAlex.TypeCrossref": "journal-article",
            "OpenAlex.Genre": "journal-article",
            "OpenAlex.Id": f"https://openalex.org/W{pmid}",
            "OpenAlex.Venue": "Test Venue",
            "OpenAlex.MeshDescriptors": "",
            "OpenAlex.MeshQualifiers": "",
            "OpenAlex.Error": "",
        }

    monkeypatch.setattr(
        "library.integration.openalex_crossref_library.fetch_openalex",
        _openalex,
    )

    def _crossref(
        session: requests.Session,
        doi: str,
        cfg: object,
        *,
        limiter: object | None = None,
        retry_cfg: object | None = None,
    ) -> dict[str, str]:
        assert isinstance(session, _StubSession)
        assert session.service == "crossref"
        return {
            "crossref.Type": "journal-article",
            "crossref.Subtype": "full-length",
            "crossref.Title": f"{doi} title",
            "crossref.Subtitle": "",
            "crossref.Subject": "testing",
            "crossref.Error": "",
        }

    monkeypatch.setattr(
        "library.integration.openalex_crossref_library.fetch_crossref",
        _crossref,
    )

    return calls


def test_document_pipeline__uses_tls_sessions(
    _patch_document_dependencies: dict[str, list[_StubSession]],
) -> None:
    pipeline = DocumentPipeline(Config())
    frame = pipeline.fetch_pubmed_records(["123456"])
    assert isinstance(frame, pd.DataFrame)
    assert not frame.empty
    assert _patch_document_dependencies["openalex"], "openalex_session was not used"
    assert _patch_document_dependencies["crossref"], "crossref_session was not used"
    assert all(session.closed for session in _patch_document_dependencies["openalex"])
    assert all(session.closed for session in _patch_document_dependencies["crossref"])


def test_cli_pubmed_mode__uses_tls_sessions(
    tmp_path: Path,
    _patch_document_dependencies: dict[str, list[_StubSession]],
) -> None:
    repo_root = Path(__file__).resolve().parents[2]
    input_path = tmp_path / "pmids.csv"
    with input_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["PMID"])
        writer.writerow(["123456"])

    output_path = tmp_path / "output.csv"
    args = [
        "--mode",
        "pubmed",
        "--config",
        str(repo_root / "config" / "config.yaml"),
        "--input",
        str(input_path),
        "--final-out",
        str(output_path),
        "--output-dir",
        str(tmp_path),
        "--limit",
        "1",
    ]

    from scripts import get_document_data as cli_module

    exit_code = cli_module.main(args)
    assert exit_code == 0
    assert output_path.exists()
    assert _patch_document_dependencies["openalex"], "openalex_session was not used"
    assert _patch_document_dependencies["crossref"], "crossref_session was not used"

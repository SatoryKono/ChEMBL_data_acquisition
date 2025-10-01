import threading
from collections import Counter
from itertools import count

import pandas as pd
import pytest
import requests

from library.config import Config
import scripts.get_document_data as document_script

from library.document_pipeline import (
    DOCUMENT_SCHEMA_COLUMNS,
    build_dataframe,
    build_quality_report,
    merge_metadata,
    merge_with_chembl,
)
from library.config import CrossRefCfg, OpenAlexCfg, PubMedCfg, SemanticScholarCfg
import scripts.get_document_data as gdd


def test_document_schema_columns_are_unique() -> None:
    """Schema column declaration should not contain duplicates."""

    assert len(DOCUMENT_SCHEMA_COLUMNS) == len(set(DOCUMENT_SCHEMA_COLUMNS))


def test_merge_metadata_normalises_fields() -> None:
    """``merge_metadata`` normalises DOI and computes publication scores."""

    pubmed = {
        "PubMed.DOI": "DOI:10.1000/ABCDEF",
        "PubMed.PublicationType": "Review",
        "PubMed.PMID": "123",
    }
    scholar = {
        "scholar.DOI": "https://doi.org/10.1000/XYZ",
        "scholar.PublicationTypes": "Comparative Study",
    }
    openalex = {"OpenAlex.PublicationTypes": "journal-article"}
    crossref = {"crossref.Type": "journal-article"}

    merged = merge_metadata(pubmed, scholar, openalex, crossref)

    assert merged["doi_normalised"] == "10.1000/abcdef"
    assert merged["doi"] == "10.1000/abcdef"
    assert merged["scholar.DOI"] == "10.1000/xyz"
    assert (
        merged["publication_types_normalised"]
        == "comparative study; journal-article; review"
    )
    assert merged["publication_type_score_review"] > 0
    assert merged["publication_type_score_experimental"] > 0
    assert merged["publication_type_score_unknown"] == 0
    assert merged["publication_class"] == "review"


def test_build_dataframe_orders_columns() -> None:
    """``build_dataframe`` respects the schema column ordering."""

    data = [{"extra": "x", "document_chembl_id": "CHEMBL1", "doi": "10/1"}]
    df = build_dataframe(data)
    expected_head = [c for c in DOCUMENT_SCHEMA_COLUMNS if c in df.columns]
    assert list(df.columns)[: len(expected_head)] == expected_head


def test_merge_with_chembl_aligns_pubmed_ids() -> None:
    """``merge_with_chembl`` joins metadata on normalised PubMed identifiers."""

    chembl_df = pd.DataFrame(
        {
            "document_chembl_id": ["CHEMBL1"],
            "pubmed_id": [123],
            "title": ["T"],
        }
    )
    meta_df = pd.DataFrame(
        {
            "document_chembl_id": ["IGNORED"],
            "title": ["Other"],
            "PubMed.PMID": ["123"],
            "PubMed.DOI": ["10.1/abc"],
        }
    )

    merged = merge_with_chembl(chembl_df, meta_df)
    assert merged["PubMed.DOI"].iloc[0] == "10.1/abc"
    assert "document_chembl_id_x" not in merged.columns
    assert "document_chembl_id_y" not in merged.columns
    assert merged["document_chembl_id"].iloc[0] == "CHEMBL1"


def test_build_quality_report_counts() -> None:
    """Quality report exposes coverage and error statistics."""

    df = pd.DataFrame(
        {
            "doi": ["10.1/abc", ""],
            "publication_class": ["review", ""],
            "PubMed.Error": ["", "err"],
            "scholar.Error": ["", ""],
            "OpenAlex.Error": ["", ""],
            "crossref.Error": ["oops", ""],
        }
    )

    report = build_quality_report(df)
    assert report["rows_total"] == 2
    assert report["doi_coverage"] == 0.5
    assert report["publication_class_counts"]["review"] == 1
    assert report["error_counts"]["pubmed"] == 1
    assert report["error_counts"]["crossref"] == 1


def test_fetch_pubmed_records_uses_fresh_sessions_per_job(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """OpenAlex and CrossRef jobs open independent retry sessions."""

    pmids = ["100", "200"]

    session_counter = count()
    counter_lock = threading.Lock()

    def _next_session_token() -> str:
        with counter_lock:
            return f"session-{next(session_counter)}"

    class DummyLimiter:
        def __init__(self, burst: int | None = None) -> None:
            self.burst = burst if burst is not None else 4

        def acquire(self) -> None:  # pragma: no cover - simple synchronisation
            return None

    pubmed_sessions: list[str] = []
    semantic_sessions: list[str] = []
    openalex_sessions: list[str] = []
    crossref_sessions: list[str] = []

    def fake_session_with_retry(*_args: object, **_kwargs: object):
        token = _next_session_token()

        class _Context:
            def __enter__(self) -> str:  # pragma: no cover - simple context
                return token

            def __exit__(
                self, *exc: object
            ) -> None:  # pragma: no cover - simple context
                return None

        return _Context()

    def fake_get_limiter(
        _name: str, _rps: float | int | None, burst: int | None = None
    ) -> DummyLimiter:
        return DummyLimiter(burst)

    def fake_pubmed_batch(
        session: str,
        batch: list[str],
        _sleep: float,
        cfg: PubMedCfg | None = None,
    ) -> list[dict[str, str]]:
        pubmed_sessions.append(session)
        return [
            {
                "PubMed.PMID": pmid,
                "PubMed.DOI": f"10.1/{pmid}",
                "PubMed.ArticleTitle": f"Article {pmid}",
                "PubMed.PublicationType": "",
            }
            for pmid in batch
        ]

    def fake_semantic_batch(
        session: str,
        batch: list[str],
        _sleep: float,
        cfg: SemanticScholarCfg | None = None,
    ) -> list[dict[str, str]]:
        semantic_sessions.append(session)
        return [
            {
                "scholar.PMID": pmid,
                "scholar.DOI": f"10.1/{pmid}",
                "scholar.PublicationTypes": "",
                "scholar.Error": "",
            }
            for pmid in batch
        ]

    def fake_semantic_single(
        *_args: object, **_kwargs: object
    ) -> dict[str, str]:  # pragma: no cover - not exercised
        return {
            "scholar.PMID": "",
            "scholar.DOI": "",
            "scholar.PublicationTypes": "",
            "scholar.Error": "",
        }

    def fake_openalex(
        session: str,
        pmid: str,
        cfg: OpenAlexCfg,
        limiter: DummyLimiter | None = None,
    ) -> dict[str, str]:
        openalex_sessions.append(session)
        return {
            "OpenAlex.PublicationTypes": "",
            "OpenAlex.TypeCrossref": "",
            "OpenAlex.Genre": "",
            "OpenAlex.Id": pmid,
            "OpenAlex.Venue": f"Venue {pmid}",
            "OpenAlex.MeshDescriptors": "",
            "OpenAlex.MeshQualifiers": "",
            "OpenAlex.Error": "",
        }

    def fake_crossref(
        session: str,
        doi: str,
        cfg: CrossRefCfg,
        limiter: DummyLimiter | None = None,
    ) -> dict[str, str]:
        crossref_sessions.append(session)
        pmid = doi.rsplit("/", 1)[-1] if doi else ""
        return {
            "crossref.Type": "",
            "crossref.Subtype": "",
            "crossref.Title": f"Title {pmid}",
            "crossref.Subtitle": "",
            "crossref.Subject": "",
            "crossref.Error": "",
        }

    monkeypatch.setattr(gdd, "session_with_retry", fake_session_with_retry)
    monkeypatch.setattr(gdd, "get_limiter", fake_get_limiter)
    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar_batch", fake_semantic_batch)
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar", fake_semantic_single)
    monkeypatch.setattr(gdd.ocl, "fetch_openalex", fake_openalex)
    monkeypatch.setattr(gdd.ocl, "fetch_crossref", fake_crossref)

    df = gdd.fetch_pubmed_records(
        pmids,
        sleep=0.0,
        pubmed_cfg=PubMedCfg(),
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=2,
        batch_size=2,
    )

    assert df["PubMed.PMID"].tolist() == pmids
    assert df["OpenAlex.Id"].tolist() == pmids
    assert df["crossref.Title"].tolist() == [f"Title {pmid}" for pmid in pmids]

    assert pubmed_sessions
    assert semantic_sessions
    assert pubmed_sessions == semantic_sessions

    assert len(openalex_sessions) == len(pmids)
    assert len(crossref_sessions) == len(pmids)
    assert len(set(openalex_sessions)) == len(pmids)
    assert len(set(crossref_sessions)) == len(pmids)

    first_session = pubmed_sessions[0]
    assert all(session != first_session for session in openalex_sessions)
    assert all(session != first_session for session in crossref_sessions)


def test_fetch_pubmed_records_reuses_duplicate_identifiers(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """OpenAlex and CrossRef lookups reuse responses for duplicate keys."""

    pmids = ["100", "200", "100"]

    class DummyLimiter:
        def __init__(self, burst: int | None = None) -> None:
            self.burst = burst if burst is not None else 4

        def acquire(self) -> None:  # pragma: no cover - simple synchronisation
            return None

    def fake_session_with_retry(*_args: object, **_kwargs: object):
        class _Context:
            def __enter__(self) -> str:  # pragma: no cover - simple context
                return "session"

            def __exit__(
                self, *_exc: object
            ) -> None:  # pragma: no cover - simple context
                return None

        return _Context()

    def fake_get_limiter(
        _name: str, _rps: float | int | None, burst: int | None = None
    ) -> DummyLimiter:
        return DummyLimiter(burst)

    def fake_pubmed_batch(
        _session: str,
        batch: list[str],
        _sleep: float,
        cfg: PubMedCfg | None = None,
    ) -> list[dict[str, str]]:
        return [
            {
                "PubMed.PMID": pmid,
                "PubMed.DOI": f"10.1/{pmid}",
                "PubMed.ArticleTitle": f"Article {pmid}",
                "PubMed.PublicationType": "",
            }
            for pmid in batch
        ]

    def fake_semantic_batch(
        _session: str,
        batch: list[str],
        _sleep: float,
        cfg: SemanticScholarCfg | None = None,
    ) -> list[dict[str, str]]:
        return [
            {
                "scholar.PMID": pmid,
                "scholar.DOI": f"10.1/{pmid}",
                "scholar.PublicationTypes": "",
                "scholar.Error": "",
            }
            for pmid in batch
        ]

    def fake_semantic_single(*_args: object, **_kwargs: object) -> dict[str, str]:
        return {
            "scholar.PMID": "",
            "scholar.DOI": "",
            "scholar.PublicationTypes": "",
            "scholar.Error": "",
        }

    openalex_calls: list[str] = []
    crossref_calls: list[str] = []

    def fake_openalex(
        _session: str,
        pmid: str,
        cfg: OpenAlexCfg,
        limiter: DummyLimiter | None = None,
    ) -> dict[str, str]:
        openalex_calls.append(pmid)
        return {
            "OpenAlex.PublicationTypes": "",
            "OpenAlex.TypeCrossref": "",
            "OpenAlex.Genre": "",
            "OpenAlex.Id": pmid,
            "OpenAlex.Venue": f"Venue {pmid}",
            "OpenAlex.MeshDescriptors": "",
            "OpenAlex.MeshQualifiers": "",
            "OpenAlex.Error": "",
        }

    def fake_crossref(
        _session: str,
        doi: str,
        cfg: CrossRefCfg,
        limiter: DummyLimiter | None = None,
    ) -> dict[str, str]:
        crossref_calls.append(doi)
        pmid = doi.rsplit("/", 1)[-1] if doi else ""
        return {
            "crossref.Type": "",
            "crossref.Subtype": "",
            "crossref.Title": f"Title {pmid}",
            "crossref.Subtitle": "",
            "crossref.Subject": "",
            "crossref.Error": "",
        }

    info_events: list[tuple[str, dict[str, object]]] = []

    def fake_info(
        event: str, *args: object, extra: dict[str, object] | None = None, **kv: object
    ) -> None:
        payload = dict(kv)
        if extra:
            payload.update(extra)
        info_events.append((event, payload))

    monkeypatch.setattr(gdd, "session_with_retry", fake_session_with_retry)
    monkeypatch.setattr(gdd, "get_limiter", fake_get_limiter)
    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar_batch", fake_semantic_batch)
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar", fake_semantic_single)
    monkeypatch.setattr(gdd.ocl, "fetch_openalex", fake_openalex)
    monkeypatch.setattr(gdd.ocl, "fetch_crossref", fake_crossref)
    monkeypatch.setattr(gdd.logger, "info", fake_info)

    df = gdd.fetch_pubmed_records(
        pmids,
        sleep=0.0,
        pubmed_cfg=PubMedCfg(),
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=2,
        batch_size=4,
    )

    assert df["PubMed.PMID"].tolist() == pmids
    assert Counter(openalex_calls) == Counter({"100": 1, "200": 1})
    assert Counter(crossref_calls) == Counter({"10.1/100": 1, "10.1/200": 1})

    cache_logs = {
        (event_payload["service"], event_payload["hits"])
        for event, event_payload in info_events
        if event == "documents_cache_reuse"
    }
    assert ("openalex", 1) in cache_logs
    assert ("crossref", 1) in cache_logs


def test_fetch_pubmed_records_logs_compact_batch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Batch failures expose compact PMID summaries in the structured log."""

    pmids = [str(100 + index) for index in range(6)]

    class DummySession:
        def __enter__(self) -> "DummySession":  # pragma: no cover - trivial context
            return self

        def __exit__(self, *_exc: object) -> None:  # pragma: no cover - trivial context
            return None

    class DummyLimiter:
        def acquire(self) -> None:  # pragma: no cover - simple synchronisation
            return None

    warning_events: list[tuple[str, dict[str, object]]] = []

    def fake_session_with_retry(*_args: object, **_kwargs: object) -> DummySession:
        return DummySession()

    def fake_get_limiter(*_args: object, **_kwargs: object) -> DummyLimiter:
        return DummyLimiter()

    def fail_pubmed_batch(*_args: object, **_kwargs: object) -> list[dict[str, str]]:
        raise requests.RequestException("boom")

    def fake_warning(
        event: str,
        *args: object,
        extra: dict[str, object] | None = None,
        **payload: object,
    ) -> None:
        record = dict(payload)
        if extra:
            record.update(extra)
        warning_events.append((event, record))

    monkeypatch.setattr(gdd, "session_with_retry", fake_session_with_retry)
    monkeypatch.setattr(gdd, "get_limiter", fake_get_limiter)
    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fail_pubmed_batch)
    monkeypatch.setattr(gdd.logger, "warning", fake_warning)

    df = gdd.fetch_pubmed_records(
        pmids,
        sleep=0.0,
        pubmed_cfg=PubMedCfg(),
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=1,
        batch_size=len(pmids),
    )

    assert len(df) == len(pmids)
    assert warning_events
    event, payload = warning_events[0]
    assert event == "pubmed_batch_request_failed"
    assert payload["error"] == "boom"
    assert payload["pmids_count"] == len(pmids)
    assert payload["pmids_sample"] == pmids[:5] + ["..."]

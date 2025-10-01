import pandas as pd

from library.config import Config
import scripts.get_document_data as document_script

from library.document_pipeline import (
    DOCUMENT_SCHEMA_COLUMNS,
    build_dataframe,
    build_quality_report,
    merge_metadata,
    merge_with_chembl,
)


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


def test_fetch_pubmed_records_skips_empty_openalex_requests(monkeypatch) -> None:
    """``fetch_pubmed_records`` should ignore blank PMIDs when querying OpenAlex."""

    cfg = Config()

    pubmed_records = [
        {"PubMed.PMID": "123", "PubMed.DOI": "10.1/abc"},
        {"PubMed.PMID": "", "PubMed.DOI": ""},
    ]

    class DummyLimiter:
        def __init__(self, burst: int | None) -> None:
            self.burst = burst if burst is not None else 1

        def acquire(self) -> None:
            return None

    class DummySession:
        def __enter__(self) -> object:
            return object()

        def __exit__(self, exc_type, exc, tb) -> bool:
            return False

    def fake_session_with_retry(*_args, **_kwargs) -> DummySession:
        return DummySession()

    def fake_get_limiter(_name: str, _rps: int, burst: int | None) -> DummyLimiter:
        return DummyLimiter(burst)

    def fake_fetch_pubmed_batch(*_args, **_kwargs) -> list[dict[str, str]]:
        return pubmed_records

    monkeypatch.setattr(document_script, "session_with_retry", fake_session_with_retry)
    monkeypatch.setattr(document_script, "get_limiter", fake_get_limiter)
    monkeypatch.setattr(document_script.pl, "fetch_pubmed_batch", fake_fetch_pubmed_batch)
    monkeypatch.setattr(
        document_script.ssl,
        "fetch_semantic_scholar_batch",
        lambda *args, **kwargs: [],
    )
    monkeypatch.setattr(
        document_script.ssl,
        "fetch_semantic_scholar",
        lambda *args, **kwargs: {},
    )

    openalex_calls: list[str] = []

    def fake_fetch_openalex(
        _session: object,
        pmid: str,
        _cfg: object,
        _limiter: object,
    ) -> dict[str, str]:
        openalex_calls.append(pmid)
        return {"OpenAlex.PMID": pmid}

    monkeypatch.setattr(document_script.ocl, "fetch_openalex", fake_fetch_openalex)
    monkeypatch.setattr(document_script.ocl, "fetch_crossref", lambda *args, **kwargs: {})

    document_script.fetch_pubmed_records(
        ["123", ""],
        sleep=0.0,
        semantic_scholar_cfg=cfg.semantic_scholar,
        openalex_cfg=cfg.openalex,
        crossref_cfg=cfg.crossref,
        pubmed_cfg=cfg.pubmed,
        max_workers=1,
        batch_size=2,
    )

    assert openalex_calls == ["123"]

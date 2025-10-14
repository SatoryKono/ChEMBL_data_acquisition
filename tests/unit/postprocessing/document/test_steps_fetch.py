from __future__ import annotations

import logging
import pandas as pd
import pytest

from library.config import Config
from library.postprocessing.document import steps


@pytest.mark.unit
@pytest.mark.parametrize(
    "raw, expected",
    [
        (" DOI:10.1000/TEST ", "10.1000/test"),
        ("https://doi.org/10.1000/ABC", "10.1000/abc"),
        (None, ""),
        ("", ""),
    ],
)
def test_clean_doi_value__normalizes_inputs(raw: object, expected: str) -> None:
    assert steps.clean_doi_value(raw) == expected


@pytest.mark.unit
def test_normalize_document_frame__coerces_columns() -> None:
    frame = pd.DataFrame(
        {
            "document_chembl_id": ["  CHEMBL123  "],
            "title": [" Example Title "],
            "doi": ["DOI:10.1111/XYZ"],
            "pubmed_id": [" 1234567 "],
            "doc_type": [" Journal  "],
            "journal": ["  Journal of Testing  "],
            "year": ["2021"],
        }
    )

    result = steps.normalize_document_frame(frame)

    assert result.loc[0, "document_chembl_id"] == "CHEMBL123"
    assert result.loc[0, "title"] == "Example Title"
    assert result.loc[0, "doi"] == "10.1111/xyz"
    assert result.loc[0, "doi_key"] == "10.1111/xyz"
    assert result.loc[0, "pubmed_id"] == "1234567"
    assert result.loc[0, "journal"] == "Journal of Testing"
    assert result.loc[0, "year"] == "2021"


@pytest.mark.unit
def test_merge_document_metadata__fills_missing_fields() -> None:
    base = pd.DataFrame(
        {
            "document_chembl_id": ["CHEMBL1"],
            "title": [pd.NA],
            "doi": ["10.2222/example"],
            "doi_key": ["10.2222/example"],
            "doc_type": [pd.NA],
            "journal": ["Journal"],
            "pubmed_id": ["123"],
            "year": ["2020"],
        }
    )
    crossref = pd.DataFrame(
        {
            "doi_key": ["10.2222/example"],
            "crossref_title": ["CrossRef Title"],
            "crossref_doc_type": ["journal-article"],
            "crossref_subject": [pd.NA],
            "crossref_error": [pd.NA],
        }
    )
    openalex = pd.DataFrame(
        {
            "doi_key": ["10.2222/example"],
            "openalex_title": ["OpenAlex Title"],
            "openalex_doc_type": ["article"],
            "openalex_type_crossref": ["journal-article"],
            "openalex_publication_year": ["2020"],
            "openalex_error": [pd.NA],
        }
    )

    merged = steps.merge_document_metadata(base, crossref_df=crossref, openalex_df=openalex)

    assert merged.loc[0, "title"] == "CrossRef Title"
    assert merged.loc[0, "doc_type"] == "journal-article"
    assert merged.loc[0, "openalex_type_crossref"] == "journal-article"


class _StubConfig(Config):
    """Lightweight configuration with deterministic defaults."""

    def __init__(self) -> None:
        super().__init__()


@pytest.mark.unit
def test_fetch_crossref_metadata__uses_fetcher_stub() -> None:
    cfg = _StubConfig()

    def _fetcher(_: str) -> tuple[dict[str, object], str]:
        return (
            {
                "message": {
                    "title": ["CrossRef Title"],
                    "type": "journal-article",
                    "subject": ["Biology", "Biology"],
                }
            },
            "",
        )

    frame = steps.fetch_crossref_metadata(
        ["10.1000/Stub"],
        config=cfg,
        fetcher=_fetcher,
    )

    assert set(frame.columns) == {
        "doi_key",
        "crossref_title",
        "crossref_doc_type",
        "crossref_subject",
        "crossref_error",
    }
    assert frame.loc[0, "crossref_title"] == "CrossRef Title"
    assert frame.loc[0, "crossref_doc_type"] == "journal-article"
    assert frame.loc[0, "crossref_subject"] == "Biology"


@pytest.mark.unit
def test_fetch_openalex_metadata__extracts_fields() -> None:
    cfg = _StubConfig()

    def _fetcher(_: str) -> tuple[dict[str, object], str]:
        return (
            {
                "doi": "https://doi.org/10.2000/test",
                "display_name": "OpenAlex Example",
                "type": "article",
                "type_crossref": "journal-article",
                "publication_year": 2023,
            },
            "",
        )

    frame = steps.fetch_openalex_metadata(
        ["123456"],
        config=cfg,
        fetcher=_fetcher,
    )

    assert set(frame.columns) == {
        "doi_key",
        "openalex_title",
        "openalex_doc_type",
        "openalex_type_crossref",
        "openalex_publication_year",
        "openalex_error",
    }
    assert frame.loc[0, "doi_key"] == "10.2000/test"
    assert frame.loc[0, "openalex_title"] == "OpenAlex Example"
    assert frame.loc[0, "openalex_doc_type"] == "article"
    assert frame.loc[0, "openalex_publication_year"] == "2023"


@pytest.mark.unit
def test_fetch_normalize_document__integrates_steps(monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = _StubConfig()
    logger = logging.getLogger("document-test")

    chembl_frame = pd.DataFrame(
        {
            "document_chembl_id": ["CHEMBL1"],
            "title": [pd.NA],
            "doi": ["10.5000/example"],
            "pubmed_id": ["987654"],
            "doc_type": [pd.NA],
            "journal": ["Journal"],
            "year": ["2022"],
        }
    )
    crossref_frame = pd.DataFrame(
        {
            "doi_key": ["10.5000/example"],
            "crossref_title": ["CrossRef Title"],
            "crossref_doc_type": ["journal"],
            "crossref_subject": [pd.NA],
            "crossref_error": [pd.NA],
        }
    )
    openalex_frame = pd.DataFrame(
        {
            "doi_key": ["10.5000/example"],
            "openalex_title": [pd.NA],
            "openalex_doc_type": ["review"],
            "openalex_type_crossref": [pd.NA],
            "openalex_publication_year": [pd.NA],
            "openalex_error": [pd.NA],
        }
    )

    monkeypatch.setattr(steps, "fetch_chembl_documents", lambda *_, **__: chembl_frame)
    monkeypatch.setattr(steps, "fetch_crossref_metadata", lambda *_, **__: crossref_frame)
    monkeypatch.setattr(steps, "fetch_openalex_metadata", lambda *_, **__: openalex_frame)

    result = steps.fetch_normalize_document(limit=1, config=cfg, logger=logger)

    assert isinstance(result, steps.DocumentFetchResult)
    assert result.data.loc[0, "title"] == "CrossRef Title"
    assert result.data.loc[0, "doc_type"] == "journal"
    assert not result.quality_report.empty
    assert isinstance(result.correlation_report, pd.DataFrame)

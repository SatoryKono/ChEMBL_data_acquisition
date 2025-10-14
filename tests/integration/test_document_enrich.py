"""Integration tests for document enrichment combining external metadata."""

from __future__ import annotations

from typing import Any

import pandas as pd
import pytest

from library.config import Config
from library.integration import openalex_crossref_library as ocl
from library.pipelines.document.pipeline import merge_metadata
from library.postprocessing.documents.steps import harmonize_document_titles


@pytest.mark.integration
def test_document_enrichment_combines_crossref_and_openalex(monkeypatch: pytest.MonkeyPatch) -> None:
    """OpenAlex and CrossRef payloads should boost coverage of document metadata."""

    calls: dict[str, int] = {"openalex": 0, "crossref": 0}

    def _fake_openalex(
        session: Any,
        pmid: str,
        *,
        cfg: Any,
        limiter: Any = None,
        retry_cfg: Any = None,
    ) -> tuple[dict[str, Any], None]:
        calls["openalex"] += 1
        return (
            {
                "id": f"https://openalex.org/W{pmid}",
                "type": "journal-article",
                "type_crossref": "journal-article",
                "genre": "journal-article",
                "host_venue": {"display_name": "OpenAlex Journal"},
                "mesh": [],
            },
            None,
        )

    def _fake_crossref(
        session: Any,
        doi: str,
        *,
        cfg: Any,
        limiter: Any = None,
        retry_cfg: Any = None,
    ) -> tuple[dict[str, Any], None]:
        calls["crossref"] += 1
        message = {
            "title": [f"CrossRef Title {doi}"],
            "type": "journal-article",
            "subtype": "article",
            "subject": ["Biology"],
        }
        return ({"message": message}, None)

    monkeypatch.setattr(ocl.openalex_client, "fetch_openalex", _fake_openalex)
    monkeypatch.setattr(ocl.crossref_client, "fetch_crossref", _fake_crossref)

    cfg = Config()
    openalex_cfg = cfg.openalex
    crossref_cfg = cfg.crossref

    sample_size = 20
    missing_index = 3
    records: list[dict[str, Any]] = []

    for idx in range(sample_size):
        pmid = f"{10_000 + idx}"
        doi = f"10.1234/example{idx}"
        chembl_title = "" if idx == missing_index else f"ChEMBL Title {idx}"
        pubmed_title = "" if idx == missing_index else f"PubMed Title {idx}"
        chembl_record = {
            "document_chembl_id": f"CHEMBL_DOC_{idx:04d}",
            "chembl.title": chembl_title,
        }
        pubmed_record = {
            "PubMed.PMID": pmid,
            "PubMed.DOI": doi,
            "PubMed.ArticleTitle": pubmed_title,
            "PubMed.PublicationType": "Journal Article" if idx % 4 else "Review",
        }
        scholar_record = {"scholar.PublicationTypes": "Journal Article"}

        openalex_payload = ocl.fetch_openalex(None, pmid, openalex_cfg)
        crossref_payload = ocl.fetch_crossref(None, doi, crossref_cfg)

        merged = merge_metadata(
            chembl_record,
            pubmed_record,
            scholar_record,
            openalex_payload,
            crossref_payload,
        )
        merged["document_chembl_id"] = chembl_record["document_chembl_id"]
        merged["chembl.title"] = chembl_title
        merged["PubMed.PMID"] = pmid
        merged["crossref.Title"] = crossref_payload.get("crossref.Title", "")
        merged["doc_type"] = merged.get("publication_class", "unknown")
        records.append(merged)

    frame = pd.DataFrame.from_records(records)
    harmonized = harmonize_document_titles(frame)

    title_series = harmonized["title"].astype("string")
    title_coverage = title_series.fillna("").str.strip().ne("").mean()
    doc_type_series = harmonized["doc_type"].astype("string")
    doc_type_coverage = doc_type_series.fillna("").str.strip().ne("").mean()

    assert title_coverage >= 0.95
    assert doc_type_coverage >= 0.95

    missing_row = harmonized.loc[harmonized["document_chembl_id"] == "CHEMBL_DOC_0003"].iloc[0]
    assert missing_row["title"] == missing_row["crossref.Title"]

    assert calls["openalex"] == sample_size
    assert calls["crossref"] == sample_size

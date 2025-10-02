"""Tests for :mod:`library.postprocessing.document`."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from library.config import IoCfg
from library.postprocessing import document as document_export_postprocessing


def _sample_dataframe() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "document_chembl_id": ["DOC1", "DOC2"],
            "title": ["Title A", ""],
            "PubMed.ArticleTitle": ["", "PubMed Title"],
            "abstract": ["", "Abstract B"],
            "PubMed.Abstract": ["Abstract A", ""],
            "doi_normalised": ["10.1000/foo", ""],
            "PubMed.DOI": ["", "10.1000/bar"],
            "pubmed_id": ["100", ""],
            "PubMed.PMID": ["", "200"],
            "year": ["2021", ""],
            "PubMed.YearCompleted": ["", "2018"],
            "journal": ["Journal A", ""],
            "PubMed.JournalTitle": ["", "Journal B"],
            "publication_class": ["experimental", "review"],
            "PubMed.is_review": [False, True],
            "OpenAlex.is_review": [False, False],
            "scholar.is_review": [False, False],
            "PubMed.MeSH_Descriptors": ["Term1; Term2", ""],
            "OpenAlex.MeshDescriptors": ["", "Term3"],
            "PubMed.MeSH_Qualifiers": ["", "Qual1"],
            "fetch_status": ["ok", "error"],
            "error_source": ["", "pubmed"],
            "PubMed.Error": ["", "Timeout"],
            "scholar.Error": ["", ""],
            "OpenAlex.Error": ["", ""],
            "crossref.Error": ["", ""],
            "date_code": ["2021-01-01", "2018-05-01"],
            "Index": ["0000", "0001"],
            "authors": ["A; B", "C"],
            "source": ["chembl", "chembl"],
            "pipeline_version": ["1.0", "1.0"],
            "timestamp_utc": ["2024-01-01T00:00:00Z", "2024-01-01T00:00:00Z"],
        }
    )


def test_preprocess_document_export_derives_fields() -> None:
    """The derived projection coalesces fields and computes coverage flags."""

    df = _sample_dataframe()
    result = document_export_postprocessing.preprocess_document_export(df)

    assert result["preferred_title"].tolist() == ["Title A", "PubMed Title"]
    assert result["preferred_abstract"].tolist() == ["Abstract A", "Abstract B"]
    assert result["preferred_doi"].tolist() == ["10.1000/foo", "10.1000/bar"]
    assert result["primary_pubmed_id"].tolist() == ["100", "200"]
    assert result["preferred_journal"].tolist() == ["Journal A", "Journal B"]
    assert result["metadata_sources"].tolist() == [
        "chembl; pubmed",
        "chembl; pubmed; openalex",
    ]
    assert result["metadata_source_count"].tolist() == [2, 3]
    assert result["has_pubmed"].tolist() == [True, True]
    assert result["has_openalex"].tolist() == [False, True]
    assert result["has_semantic_scholar"].tolist() == [False, False]
    assert result["error_sources"].tolist() == ["", "pubmed"]
    assert result["has_error"].tolist() == [False, True]
    assert result["mesh_terms"].tolist() == ["Term1; Term2", "Qual1; Term3"]
    assert result["publication_year"].tolist() == [2021, 2018]
    assert result["is_review"].tolist() == [False, True]


def test_postprocess_export_file_writes_csv(tmp_path: Path) -> None:
    """``postprocess_export_file`` writes the derived projection to disk."""

    df = _sample_dataframe()
    cfg = IoCfg(csv_sep=",", csv_encoding="utf-8")
    input_path = tmp_path / "output.documents_20240101.csv"
    df.to_csv(input_path, index=False, sep=cfg.csv_sep, encoding=cfg.csv_encoding)

    output_path = document_export_postprocessing.postprocess_export_file(
        input_path,
        cfg=cfg,
    )

    assert output_path.name == f"preprocessed_{input_path.name}"
    processed = pd.read_csv(output_path, sep=cfg.csv_sep, encoding=cfg.csv_encoding)
    assert "metadata_sources" in processed.columns
    assert processed["metadata_sources"].iloc[1] == "chembl; pubmed; openalex"
